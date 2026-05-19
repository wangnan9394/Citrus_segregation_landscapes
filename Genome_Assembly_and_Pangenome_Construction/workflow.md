# *Citrus hindsii* Pangenome Assembly and Analysis Workflow


## Step 1 — Haplotype-Resolved Genome Assembly

We performed haplotype-resolved *de novo* assembly for PN and DB using PacBio HiFi long reads combined with Hi-C chromatin interaction data via hifiasm v0.19.6. This integrative approach simultaneously phases the two haplotypes of each parental line based on Hi-C contact signals, producing per-haplotype contig-level assemblies in GFA format. Default parameters were applied as recommended for diploid genomes with high heterozygosity.

```bash
for example PN
/home/wangnan/Software/hifiasm-0.19.5/hifiasm -t 20 -o PNasm --h1 /home/wangnan/8.kumquat/02.HIFI_HIC_data/HIC/PN02_R1.fq.gz --h2 /home/wangnan/8.kumquat/02.HIFI_HIC_data/HIC/PN02_R2.fq.gz  /home/wangnan/8.kumquat/02.HIFI_HIC_data/HIFI/PN02.CCS.fq.gz
```

**Output**: PNasm.hic.hap1.p_ctg.gfa, PNasm.hic.hap2.p_ctg.gfa, DBasm.hic.hap1.p_ctg.gfa, DBasm.hic.hap2.p_ctg.gfa

---

## Step 2 — GFA to FASTA Conversion

The hifiasm assembly graph outputs are in GFA format, which must be converted to standard FASTA for downstream scaffolding and annotation tools. gfatools extracts the contig sequences embedded in each haplotype GFA and writes them as plain FASTA records.

```bash
for i in 1 2; do gfatools gfa2fa ${PREFIX}.hic.hap${i}.p_ctg.gfa > hap${i}.p_ctg.fa; done
```

**Output**: `PNasm.hic.hap1.p_ctg.fa`, `PNasm.hic.hap2.p_ctg.fa` (and equivalently for DB)

---

## Step 3 — Organelle Contig Filtering

Assembled contig sets routinely include mitochondrial and chloroplast sequences, which must be removed before nuclear genome scaffolding to prevent spurious Hi-C contact signals and inflated assembly statistics. filter_organelle.py aligns each contig against reference organellar genomes and flags those with high similarity for exclusion, retaining only nuclear-origin contigs.

```bash
python contigs.filter.py -f ${i} -c CP.fa -m MITO.fa -t 20 -o ${i}_result.tx
```

```filter_organelle.py
#!/usr/bin/env python3
"""
优化版基因组contig筛选脚本
功能：从组装结果中筛选高质量的核基因组contigs，过滤细胞器污染
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import shutil
from Bio import SeqIO
import gzip
import logging
from datetime import datetime

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(f'genome_filter_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)


class GenomeFilter:
    def __init__(self, args):
        """初始化参数"""
        self.fasta = Path(args.fasta)
        self.cp_ref = Path(args.req_cp)
        self.mt_ref = Path(args.req_mt)
        self.threads = int(args.threads)
        self.output = Path(args.output)
        self.work_dir = Path(tempfile.mkdtemp(prefix="genome_filter_"))
        
        # 筛选阈值
        self.min_length = 100000
        self.max_gc = 50.0
        self.max_cp_coverage = 0.01
        self.max_mt_coverage = 0.01
        self.max_unique_kmer = 0.7
        
        # 检查输入文件
        self._validate_inputs()
        
        logger.info(f"工作目录: {self.work_dir}")
        logger.info(f"输入文件: {self.fasta}")
        logger.info(f"使用线程数: {self.threads}")
    
    def _validate_inputs(self):
        """验证输入文件是否存在"""
        required_files = [
            (self.fasta, "输入FASTA文件"),
            (self.cp_ref, "叶绿体参考基因组"),
            (self.mt_ref, "线粒体参考基因组")
        ]
        
        for filepath, desc in required_files:
            if not filepath.exists():
                logger.error(f"{desc}不存在: {filepath}")
                raise FileNotFoundError(f"{desc}不存在: {filepath}")
            logger.info(f"找到{desc}: {filepath}")
    
    def _parse_fasta(self) -> Dict[str, str]:
        """解析FASTA文件，返回序列字典"""
        sequences = {}
        try:
            with gzip.open(self.fasta, 'rt') if str(self.fasta).endswith('.gz') else open(self.fasta, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    sequences[record.id] = str(record.seq)
            logger.info(f"成功解析 {len(sequences)} 条序列")
            return sequences
        except Exception as e:
            logger.error(f"解析FASTA文件失败: {e}")
            raise
    
    def calculate_basic_stats(self, sequences: Dict[str, str]) -> Tuple[Dict, Dict]:
        """计算序列长度和GC含量"""
        lengths = {}
        gc_contents = {}
        
        for seq_id, seq in sequences.items():
            lengths[seq_id] = len(seq)
            gc_count = seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')
            gc_contents[seq_id] = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0
        
        logger.info("基本统计计算完成")
        return lengths, gc_contents
    
    def _run_minimap2(self, query: Path, target: Path, output_paf: Path, target_name: str) -> None:
        """运行minimap2进行比对"""
        cmd = [
            'minimap2', '-t', str(self.threads),
            '-x', 'asm5',
            str(target), str(query),
            '-o', str(output_paf)
        ]
        
        logger.info(f"运行minimap2比对到{target_name}...")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            if result.stderr:
                logger.debug(f"minimap2 {target_name} stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"minimap2运行失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            raise
    
    def calculate_coverage(self, sequences: Dict[str, str], lengths: Dict[str, int], 
                          target_file: Path, target_name: str) -> Dict[str, float]:
        """计算与参考基因组的覆盖度"""
        # 创建临时查询文件
        query_file = self.work_dir / f"query_{target_name}.fa"
        with open(query_file, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
        
        # 运行minimap2
        paf_file = self.work_dir / f"alignment_{target_name}.paf"
        self._run_minimap2(query_file, target_file, paf_file, target_name)
        
        # 解析PAF文件计算覆盖度
        coverage_dict = {}
        coverage_counts = {seq_id: 0 for seq_id in sequences.keys()}
        
        try:
            with open(paf_file, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 11:
                        seq_id = fields[0]
                        if seq_id in coverage_counts:
                            start = int(fields[2])
                            end = int(fields[3])
                            coverage_counts[seq_id] += (end - start)
            
            # 计算覆盖度比例
            for seq_id, covered_bases in coverage_counts.items():
                if lengths[seq_id] > 0:
                    coverage_dict[seq_id] = covered_bases / lengths[seq_id]
                else:
                    coverage_dict[seq_id] = 0.0
            
            logger.info(f"{target_name}覆盖度计算完成")
            return coverage_dict
            
        except Exception as e:
            logger.error(f"解析{target_name}比对结果失败: {e}")
            raise
    
    def _run_meryl_single(self, seq_id: str, seq: str) -> Tuple[str, float]:
        """为单条序列运行meryl计算唯一k-mer"""
        # 创建临时序列文件
        seq_file = self.work_dir / f"{seq_id}.fa"
        meryl_db = self.work_dir / f"{seq_id}.meryl"
        
        with open(seq_file, 'w') as f:
            f.write(f">{seq_id}\n{seq}\n")
        
        try:
            # 运行meryl count
            count_cmd = [
                'meryl', f'threads={min(4, self.threads)}',
                'count', 'compress', 'k=41',
                str(seq_file), 'output', str(meryl_db)
            ]
            
            subprocess.run(count_cmd, capture_output=True, check=True)
            
            # 获取统计信息
            stats_cmd = ['meryl', 'statistics', str(meryl_db)]
            result = subprocess.run(stats_cmd, capture_output=True, text=True, check=True)
            
            # 解析结果获取唯一k-mer数
            unique_kmers = 0
            for line in result.stdout.split('\n'):
                if 'distinct' in line and 'non-redundant' in line and 'missing' not in line:
                    parts = line.replace(' ', '').split('(')[0]
                    unique_kmers = int(parts.split('t')[-1])
                    break
            
            # 计算唯一k-mer比例
            unique_ratio = 1 - (unique_kmers / len(seq)) if len(seq) > 0 else 0
            
            # 清理临时文件
            seq_file.unlink(missing_ok=True)
            shutil.rmtree(meryl_db, ignore_errors=True)
            
            return seq_id, unique_ratio
            
        except Exception as e:
            logger.error(f"计算序列 {seq_id} 的唯一k-mer时出错: {e}")
            # 清理临时文件
            seq_file.unlink(missing_ok=True)
            shutil.rmtree(meryl_db, ignore_ok=True)
            return seq_id, 1.0  # 出错时返回最大值
    
    def calculate_unique_kmers(self, sequences: Dict[str, str], lengths: Dict[str, int]) -> Dict[str, float]:
        """并行计算所有序列的唯一k-mer比例"""
        logger.info("开始计算唯一k-mer比例...")
        
        unique_ratios = {}
        
        # 使用进程池并行处理
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = {
                executor.submit(self._run_meryl_single, seq_id, seq): seq_id
                for seq_id, seq in sequences.items()
            }
            
            completed = 0
            total = len(futures)
            for future in as_completed(futures):
                seq_id, ratio = future.result()
                unique_ratios[seq_id] = ratio
                completed += 1
                if completed % 10 == 0 or completed == total:
                    logger.info(f"进度: {completed}/{total} 条序列处理完成")
        
        logger.info("唯一k-mer比例计算完成")
        return unique_ratios
    
    def filter_sequences(self, sequences: Dict[str, str], lengths: Dict[str, int], 
                        gc_contents: Dict[str, float], cp_coverage: Dict[str, float],
                        mt_coverage: Dict[str, float], unique_ratios: Dict[str, float]) -> List[str]:
        """根据条件筛选序列"""
        selected = []
        
        for seq_id in sequences.keys():
            length_ok = lengths[seq_id] > self.min_length
            gc_ok = gc_contents[seq_id] < self.max_gc
            cp_ok = cp_coverage[seq_id] < self.max_cp_coverage
            mt_ok = mt_coverage[seq_id] < self.max_mt_coverage
            unique_ok = unique_ratios[seq_id] < self.max_unique_kmer
            
            if all([length_ok, gc_ok, cp_ok, mt_ok, unique_ok]):
                selected.append(seq_id)
        
        logger.info(f"筛选结果: 共 {len(sequences)} 条序列，{len(selected)} 条通过筛选")
        return selected
    
    def write_output(self, sequences: Dict[str, str], lengths: Dict[str, int], 
                    gc_contents: Dict[str, float], cp_coverage: Dict[str, float],
                    mt_coverage: Dict[str, float], unique_ratios: Dict[str, float],
                    selected_ids: List[str]) -> None:
        """写入输出文件"""
        # 写入统计文件
        with open(self.output, 'w') as f:
            # 写入表头
            f.write("Seq_ID\tLength\tGC(%)\tCP_Coverage\tMT_Coverage\tUnique_kmer_ratio\tSelected\n")
            
            for seq_id in sequences.keys():
                selected = "Yes" if seq_id in selected_ids else "No"
                f.write(f"{seq_id}\t{lengths[seq_id]}\t{gc_contents[seq_id]:.2f}\t"
                       f"{cp_coverage[seq_id]:.6f}\t{mt_coverage[seq_id]:.6f}\t"
                       f"{unique_ratios[seq_id]:.6f}\t{selected}\n")
        
        # 写入筛选后的FASTA文件
        clean_fasta = self.output.parent / f"{self.output.stem}.CLEAN.fa"
        with open(clean_fasta, 'w') as f:
            for seq_id in selected_ids:
                f.write(f">{seq_id}\n{sequences[seq_id]}\n")
        
        logger.info(f"统计结果已写入: {self.output}")
        logger.info(f"筛选后的FASTA已写入: {clean_fasta}")
        
        # 生成筛选统计摘要
        summary_file = self.output.parent / f"{self.output.stem}.summary.txt"
        with open(summary_file, 'w') as f:
            f.write("=== 基因组过滤结果摘要 ===\n")
            f.write(f"输入序列总数: {len(sequences)}\n")
            f.write(f"筛选后序列数: {len(selected_ids)}\n")
            f.write(f"保留比例: {len(selected_ids)/len(sequences)*100:.2f}%\n")
            f.write(f"总碱基数: {sum(lengths.values()):,}\n")
            f.write(f"筛选后碱基数: {sum(lengths[s] for s in selected_ids):,}\n")
            f.write("\n=== 筛选条件 ===\n")
            f.write(f"最小长度: {self.min_length:,} bp\n")
            f.write(f"最大GC含量: {self.max_gc}%\n")
            f.write(f"最大叶绿体覆盖度: {self.max_cp_coverage*100:.1f}%\n")
            f.write(f"最大线粒体覆盖度: {self.max_mt_coverage*100:.1f}%\n")
            f.write(f"最大唯一k-mer比例: {self.max_unique_kmer*100:.1f}%\n")
        
        logger.info(f"统计摘要已写入: {summary_file}")
    
    def run(self):
        """主运行函数"""
        logger.info("开始基因组过滤分析...")
        
        try:
            # 1. 解析FASTA文件
            sequences = self._parse_fasta()
            
            # 2. 计算基本统计信息
            lengths, gc_contents = self.calculate_basic_stats(sequences)
            
            # 3. 计算叶绿体覆盖度
            cp_coverage = self.calculate_coverage(sequences, lengths, self.cp_ref, "CP")
            
            # 4. 计算线粒体覆盖度
            mt_coverage = self.calculate_coverage(sequences, lengths, self.mt_ref, "MT")
            
            # 5. 计算唯一k-mer比例
            unique_ratios = self.calculate_unique_kmers(sequences, lengths)
            
            # 6. 筛选序列
            selected_ids = self.filter_sequences(
                sequences, lengths, gc_contents, 
                cp_coverage, mt_coverage, unique_ratios
            )
            
            # 7. 写入输出
            self.write_output(
                sequences, lengths, gc_contents,
                cp_coverage, mt_coverage, unique_ratios,
                selected_ids
            )
            
            logger.info("分析完成！")
            
        except Exception as e:
            logger.error(f"分析过程中发生错误: {e}", exc_info=True)
            raise
        
        finally:
            # 清理工作目录
            if self.work_dir.exists():
                shutil.rmtree(self.work_dir)
                logger.info(f"已清理临时工作目录: {self.work_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='基因组contig筛选工具 - 从组装结果中筛选高质量的核基因组contigs，过滤细胞器污染',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python genome_filter.py -f contigs.fa -c cp.fasta -m mt.fasta -t 16 -o results.tsv
        
输出文件:
  results.tsv - 包含所有序列的详细统计信息
  results.CLEAN.fa - 筛选后的序列
  results.summary.txt - 分析摘要
        
筛选条件:
  长度 > 100,000 bp
  GC含量 < 50%
  叶绿体覆盖度 < 1%
  线粒体覆盖度 < 1%
  唯一k-mer比例 < 70%
        """
    )
    
    parser.add_argument('-f', '--fasta', required=True, 
                       help='输入FASTA文件 (支持.gz压缩格式)')
    parser.add_argument('-c', '--req_cp', required=True, 
                       help='叶绿体参考基因组文件')
    parser.add_argument('-m', '--req_mt', required=True, 
                       help='线粒体参考基因组文件')
    parser.add_argument('-t', '--threads', required=True, type=int,
                       help='使用的线程数')
    parser.add_argument('-o', '--output', required=True,
                       help='输出统计文件路径')
    
    # 可选参数：调整筛选阈值
    parser.add_argument('--min_length', type=int, default=100000,
                       help='最小序列长度 (默认: 100000)')
    parser.add_argument('--max_gc', type=float, default=50.0,
                       help='最大GC含量百分比 (默认: 50.0)')
    parser.add_argument('--max_cp_cov', type=float, default=0.01,
                       help='最大叶绿体覆盖度 (默认: 0.01)')
    parser.add_argument('--max_mt_cov', type=float, default=0.01,
                       help='最大线粒体覆盖度 (默认: 0.01)')
    parser.add_argument('--max_unique_kmer', type=float, default=0.7,
                       help='最大唯一k-mer比例 (默认: 0.7)')
    
    args = parser.parse_args()
    
    try:
        filter = GenomeFilter(args)
        filter.run()
        
    except KeyboardInterrupt:
        logger.info("用户中断程序执行")
        sys.exit(1)
    except Exception as e:
        logger.error(f"程序执行失败: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
```

---

## Step 4 — Hi-C Scaffolding

### 4.1 Initial Per-haplotype Scaffolding (Juicer + 3D-DNA + Juicebox)

Hi-C reads carry long-range chromatin contact information that can anchor and orient contigs into chromosome-scale scaffolds. For each haplotype, Hi-C reads were aligned with Juicer v1.6 to generate contact matrices; 3D-DNA v201008 then clustered, ordered and oriented contigs into pseudo-chromosomes. Manual curation in Juicebox v2.20.0 was performed to correct mis-joins and refine chromosome boundaries based on haplotype-specific interaction heatmaps.

```bash
bash /home/wangnan/Software/juicer-main/scripts/juicer.sh -d $PWD -D /home/wangnan/Software/juicer-main/ -z ./references/*.fa -y ./restriction_sites/*_DpnII.txt -p ./restriction_sites/*.chrom.sizes -s DpnII -t 16 -g PN
samtools view -@ 16 -O SAM -F 1024 ./aligned/merged_dedup.bam | awk -v mnd=60 -f /home/wangnan/Software/juicer-main/scripts/common/sam_to_pre.awk > 60_F1024_merged_nodups.txt
mkdir 3d-dna && cd 3d-dna
bash /home/wangnan/Software/3d-dna-master/run-asm-pipeline.sh -r 2 ../references/*.fa ../60_F1024_merged_nodups.txt
mkdir review-3d-dna && cd review-3d-dna
bash /home/wangnan/Software/3d-dna-master/run-asm-pipeline-post-review.sh -r ../PN.haps.0.review.assembly ../references/*.fa ../0_q0_merged_nodups.txt
```

---

### 4.2 Cross-haplotype Integration and HapHiC Refinement

After individual scaffolding, the two chromosome-scale haplotypes of each parent (PN_hap1 + PN_hap2, or DB_hap1 + DB_hap2) were merged into a single integrated dataset. Hi-C reads were then remapped to this combined assembly using HapHiC v1.0.7, which leverages haplotype-aware contact information to correct interaction errors and further improve scaffolding accuracy. This two-step strategy produces the final haplotype-resolved genome for each parent.

```bash
# Concatenate two haplotypes of the same parent
cat PN_hap1.nuclear.FINAL.fasta PN_hap2.nuclear.FINAL.fasta > PN_integrated.fa

# HapHiC re-scaffolding for accuracy refinement
haphic pipeline PN_integrated.fa ${HIC_BAM} 2 --threads ${THREADS} --outdir haphic_PN_out
```

---

## Step 5 — Genome Completeness Assessment

Genome completeness of each haplotype assembly was evaluated with BUSCO v5.8.3 against the embryophyta_odb10 database. High BUSCO completeness (>95% complete BUSCOs) confirms sufficient assembly quality for downstream annotation and pangenome analyses.

```bash
busco -i hap${i}.final.fa -l embryophyta_odb10 -o busco_hap${i} -m genome --cpu ${THREADS}
```

---

## Step 6 — Genome Annotation

### 6.1 Transposable Element (TE) Annotation

Repetitive element annotation was performed using EDTA v2.1.0, which builds a species-specific TE library through a comprehensive pipeline covering LTR retrotransposons, TIR transposons, and other repeat classes. The resulting TE library was subsequently used for repeat masking prior to gene prediction.

```bash
EDTA.pl --genome hap${i}.final.fa --species others --step all --sensitive 1 --anno 1 --threads ${THREADS}
```

**Output**: `hap${i}.final.fa.mod.EDTA.TEanno.gff3`, `hap${i}.final.fa.mod.EDTA.TElib.fa`

---

### 6.3 RNA-seq Alignment (Transcriptome Evidence)

Previously published RNA-seq reads from multiple *C. hindsii* tissues (root, stem, leaf, flower and fruit; PRJNA735863) were quality-filtered and aligned to the repeat-masked genome using HISAT2. The resulting alignments provide transcript-level evidence to guide gene model building in PASA and EVM.

```bash
hisat2 -x hap${i}.masked.index -1 ${RNA_R1} -2 ${RNA_R2} -p ${THREADS} | samtools sort -o hap${i}.rna.bam
```

**Output**: `hap${i}.rna.bam`

---

### 6.4 Gene Model Integration (GeMoMa + PASA + AUGUSTUS + EVM)

Gene models were built by integrating multiple lines of evidence: homologous protein alignments via GeMoMa v0.98.5, transcript assemblies via PASA v2.5.3, *de novo* structure predictions via AUGUSTUS v3.5.0, and final consensus integration via EVidenceModeler v2.1.0. This multi-evidence strategy improves gene model accuracy and completeness relative to any single method.

```bash
# Step 1: Homology-based prediction with GeMoMa
GeMoMa -Xmx64G GeMoMaPipeline tseq=hap${i}.final.fa.masked a=protein_db.fa o=true AnnotationFinalizer.r=NO t=${THREADS}

# Step 2: Transcript assembly with PASA
Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g hap${i}.final.fa.masked -t transcripts.fa --ALIGNERS blat,gmap --CPU ${THREADS}

# Step 3: De novo prediction with AUGUSTUS
augustus --species=Citrus_MODEL --softmasking=1 --gff3=on hap${i}.final.fa.masked > augustus_hap${i}.gff3

# Step 4: Consensus integration with EVM
EVidenceModeler --sample_id hap${i} --genome hap${i}.final.fa.masked --gene_predictions gemoma.gff3,augustus_hap${i}.gff3 --transcript_alignments pasa.gff3 --repeats hap${i}.final.fa.out.gff --CPU ${THREADS}
```

---

## Step 7 — Per-chromosome Sequence Extraction

The four phased haplotypes (PN_hap1, PN_hap2, DB_hap1, DB_hap2) are split by individual chromosome in preparation for chromosome-level pangenome graph construction. Splitting by chromosome reduces computational complexity and makes graph building tractable on standard HPC nodes.

```bash
for CHR in $(seq -w 1 9 | awk '{print "chr"$1}'); do for i in PN_hap1 PN_hap2 DB_hap1 DB_hap2; do samtools faidx ${i}.final.fa ${CHR} > pangenome/${CHR}/${i}.${CHR}.fa; done; done
```

---

## Step 8 — Pangenome Graph Construction

### 8.1 PGGB — Graph Visualization and Structural Variation

PGGB v0.7.3 builds a sequence-to-sequence alignment graph from all four haplotypes simultaneously for each chromosome. This all-against-all graph encodes complex structural variation including inversions, duplications and large indels as graph topology, enabling visualization of haplotype paths through `odgi` and detection of SVs without a linear reference bias.

```bash
for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9; do pggb -i pangenome/${CHR}/all_haps.${CHR}.fa.gz -p 90 -s 10000 -n 117 -k 47 -P 1,4,6,2,26,1 -O 0.001  -t 40 -V 'SJG' -o pangenome/${CHR}/pggb_out; done
```

---

### 8.2 Minigraph-Cactus (MC) — Variant Calling and VCF Generation

Minigraph-Cactus v3.0.0 constructs a hierarchical pangenome graph using a reference-guided progressive alignment strategy. The *C. hindsii* reference genome v2.0 serves as the backbone anchor, and the four haplotypes are progressively integrated. MC directly outputs a VCF representing all structural and sequence variants across haplotypes, enabling standard variant filtering and genotyping workflows. The `--vcf` flag activates VCF deconstruction during graph construction.

```bash
/home/wangnan/Software/cactus-bin-v2.9.8/venv-cactus-v2.9.8/bin/cactus-pangenome js haps.txt --maxCores 20 --maxMemory 800G --maxDisk 10000G --outDir Pangenome --outName KumquatPG --reference SJG --vcf --giraffe --gfa --gbz --odgi --workDir $PWD --odgi --workDir $PWD --restart
```

---
