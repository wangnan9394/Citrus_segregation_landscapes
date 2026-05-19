# *Citrus hindsii* PN × DB F1 mapping population (226 offspring)  

## Step 1 — Collinearity Analysis and Large-Scale Inversion Identification

Synteny analysis was performed among the four haplotype-resolved genomes using the GENESPACE pipeline v1.2.3, which infers genome-wide collinear blocks from annotated gene sets via OrthoFinder and MCScanX. High-confidence protein sequences extracted from each haplotype annotation served as input. Collinear patterns and structural rearrangements were visualized with LINKVIEW2; megabase-level inversions were systematically screened and their breakpoints and genomic coordinates documented for downstream characterization.

```r
# Step 1a: Extract longest-isoform protein per gene
# Step 1b: Run GENESPACE (R)
#cat DBhap1_gene.space.work.r 
library(GENESPACE)
gpar <- init_genespace(
  wd = "/home/wangnan/8.kumquat/05.assemble/05.helixer_genespace_pre/01.genespace/DBhap1",
  path2mcscanx = "/home/wangnan/Software/MCScanX",
  path2orthofinder = "/home/wangnan/Software/OrthoFinder_source/orthofinder",
  path2diamond = "/home/wangnan/ENV/miniconda3/envs/genespace/bin/diamond")
out <- run_genespace(gpar,overwrite=F)"
```
```bash
# Step 1c: Conduct GENESPACE collinear blocks
Rscript ${i}_gene.space.work.r

# Step 1d: Visualize with LINKVIEW2
LINKVIEW2 --bezier --style simple -k k.txt --chro_axis --svg_width 1800 --svg_height 200 --no_label --svg_space 0 --align left --gap_length 0.01 USED_SYS.DB2_vs_DB1.synHits.txt -o DB.plot
LINKVIEW2 --bezier --style simple -k k.txt --chro_axis --svg_width 1800 --svg_height 200 --no_label --svg_space 0 --align left --gap_length 0.01 USED_SYS.PN2_vs_PN1.synHits.txt -o PN.plot
```

---

## Step 2 — Haplotype-Specific Variant Map Construction

Raw variants from the Minigraph-Cactus pangenome VCF were filtered with two stringent criteria to generate a high-confidence parental haplotype variant atlas: (i) the variant site must be captured by all four haplotypes; (ii) markers were classified by their relationship to the reference backbone — those concordant with the reference (e.g., genotype pattern 0;1;1;1) and those discordant (e.g., 1;0;0;0). The resulting atlas contains haplotype-specific SNPs, InDels and presence-absence variations (PAVs) for each of the four haplotypes, and serves as the marker foundation for population-level haplotype tracing.

```bash
python filter.hap.specific.py
```
```python
#filter.hap.specific.py
outputDB1=open('DBhap1.specific.1.vcf','w')
outputDB2=open('DBhap2.specific.1.vcf','w')
outputPN1=open('PNhap1.specific.1.vcf','w')
outputPN2=open('PNhap2.specific.1.vcf','w')
with open("KumquatPG.snp.two_alleles.vcf","r") as fd:
    for line0 in fd:
        if line0.startswith("#"):
            outputDB1.write(line0)
            outputDB2.write(line0)
            outputPN1.write(line0)
            outputPN2.write(line0)
            continue
        line = line0.strip().replace("\n","").split("\t")
        db1 = line[9]
        db2 = line[10]
        pn1 = line[11]
        pn2 = line[12]
        if db1 == '1' and db2 == '0' and pn1 == '0' and pn2 == '0':
            outputDB1.write(line0)
        if db1 == '0' and db2 == '1' and pn1 == '0' and pn2 == '0':
            outputDB2.write(line0)
        if db1 == '0' and db2 == '0' and pn1 == '1' and pn2 == '0':
            outputPN1.write(line0)
        if db1 == '0' and db2 == '0' and pn1 == '0' and pn2 == '1':
            outputPN2.write(line0)
        



outputDB1=open('DBhap1.specific.0.vcf','w')
outputDB2=open('DBhap2.specific.0.vcf','w')
outputPN1=open('PNhap1.specific.0.vcf','w')
outputPN2=open('PNhap2.specific.0.vcf','w')
with open("KumquatPG.snp.two_alleles.vcf","r") as fd:
    for line0 in fd:
        if line0.startswith("#"):
            outputDB1.write(line0)
            outputDB2.write(line0)
            outputPN1.write(line0)
            outputPN2.write(line0)
            continue
        line = line0.strip().replace("\n","").split("\t")
        db1 = line[9]
        db2 = line[10]
        pn1 = line[11]
        pn2 = line[12]
        if db1 == '0' and db2 == '1' and pn1 == '1' and pn2 == '1':
            outputDB1.write(line0)
        if db1 == '1' and db2 == '0' and pn1 == '1' and pn2 == '1':
            outputDB2.write(line0)
        if db1 == '1' and db2 == '1' and pn1 == '0' and pn2 == '1':
            outputPN1.write(line0)
        if db1 == '1' and db2 == '1' and pn1 == '1' and pn2 == '0':
            outputPN2.write(line0)
```

---

## Step 3 — Population Genotyping with the Pangenome Graph (vg v1.70.0)

Illumina resequencing reads from 226 F1 hybrid offspring were aligned to the pangenome graph using `vg giraffe`, which maps short reads directly onto the graph index to avoid reference bias. Population-level variants were then called per individual using `vg call`. The per-sample VCFs were intersected with the parental haplotype marker atlas to extract haplotype-informative genotypes, generating a population-wide haplotype genotype matrix for downstream crossover and segregation analyses.

```bash
threads=8
ref_dir=/home/wangnan/8.kumquat/08.pangenome_based_alignment
nbt_ngs=/home/wangnan/8.kumquat/03.Fh_genetic_pops
cd /home/wangnan/8.kumquat/08.pangenome_based_alignment/1.results
vg giraffe -p -r 1 -b fast -t \$threads -Z \$ref_dir/merge_call.giraffe.gbz -d \$ref_dir/merge_call.dist -m \$ref_dir/merge_call.min -f \$nbt_ngs/${i}_1_clean.fq.gz -f \$nbt_ngs/${i}_2_clean.fq.gz > ${i}.gam
vg stats -a ${i}.gam > ${i}.stats
vg pack -x \$ref_dir/merge_call.giraffe.gbz -g ${i}.gam -o ${i}.pack -Q 5
vg call -t \$threads \$ref_dir/merge_call.giraffe.gbz -r \$ref_dir/merge_call.snarls -k ${i}.pack -s ${i} > ${i}.vcf
bgzip ${i}.vcf
bcftools merge -m all --threads $threads *ALL_VCF_gz_files* -Oz -o NEW.merged_queries.vcf.gz

python step1.基于亲本过滤vcf.py input.vcf output.vcf
python step2.转换vcf为新的hap标记.py
```
```python
#step1.基于亲本过滤vcf.py
import sys

def filter_vcf_simple(input_file: str, output_file: str) -> None:
    """
    简化版VCF过滤
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        db_idx = -1
        pn_idx = -1
        
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    try:
                        db_idx = headers.index('DB')
                        pn_idx = headers.index('PN')
                        print(f"找到DB列: 索引{db_idx}, PN列: 索引{pn_idx}")
                    except ValueError as e:
                        print(f"错误: 未找到样本列 - {e}")
                        sys.exit(1)
                continue
            
            cols = line.strip().split('\t')
            if len(cols) <= max(db_idx, pn_idx):
                continue
            
            # 提取基因型（假设GT是第一个字段）
            db_gt = cols[db_idx].split(':')[0] if ':' in cols[db_idx] else cols[db_idx]
            pn_gt = cols[pn_idx].split(':')[0] if ':' in cols[pn_idx] else cols[pn_idx]
            
            # 过滤条件
            if db_gt == '0/1' and pn_gt == '0/0':
                f_out.write(line)
            if db_gt == '1/0' and pn_gt == '0/0':
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("使用方法: python filter_vcf.py input.vcf output.vcf")
        sys.exit(1)
    
    filter_vcf_simple(sys.argv[1], sys.argv[2])
    print("过滤完成!")
```
```python
#step2.转换vcf为新的hap标记.py
output=open("DBhap1.specific.1.vcf.POS.filtered.DB.PN.vcf.HAPS.txt","w")
with open("DBhap1.specific.1.vcf.POS.filtered.DB.PN.vcf","r") as fd:
    for line0 in fd:
        if line0.startswith("#CHR"):
            output.write(line0)
        if not line0.startswith("#"):
            line = line0.strip().replace("\n","").split("\t")
            #print(line)
            early_part=line[:9]
            gt=line[10:-1]
            new_tag=[]
            for each in gt:
                if each.split(':')[0]== '0/0':
                    tag="DBhap2"
                    new_tag.append(tag)
                elif each.split(':')[0]== '1/0' or each.split(':')[0]== '0/1':
                    tag="DBhap1"
                    new_tag.append(tag)
                else:
                    tag='NA'
                    new_tag.append(tag)
            ll="\t".join(early_part)+'\t'+"\t".join(new_tag)+'\n'
            output.write(ll)
            #break
output.close()

output=open("DBhap2.specific.1.vcf.POS.filtered.DB.PN.vcf.HAPS.txt","w")
with open("DBhap2.specific.1.vcf.POS.filtered.DB.PN.vcf","r") as fd:
    for line0 in fd:
        if line0.startswith("#CHR"):
            output.write(line0)
        if not line0.startswith("#"):
            line = line0.strip().replace("\n","").split("\t")
            #print(line)
            early_part=line[:9]
            gt=line[10:-1]
            new_tag=[]
            for each in gt:
                if each.split(':')[0]== '0/0':
                    tag="DBhap1"
                    new_tag.append(tag)
                elif each.split(':')[0]== '1/0' or each.split(':')[0]== '0/1':
                    tag="DBhap2"
                    new_tag.append(tag)
                else:
                    tag='NA'
                    new_tag.append(tag)
            ll="\t".join(early_part)+'\t'+"\t".join(new_tag)+'\n'
            output.write(ll)
            #break
```

---

## Step 4 — Crossover Identification

### 4.1 Sliding-window Haplotype Purity Scoring

Recombination breakpoints across all 226 F1 offspring were identified using a sliding-window strategy applied to the haplotype genotype matrix. Windows of 20 consecutive haplotype-specific markers were scored: a window was classified as "pure" if ≥ 90% of its markers derived from a single parental haplotype (hap1 or hap2), and as "mixed" when the marker proportions deviated from this purity threshold, indicating a potential crossover region. This density-aware approach tolerates sporadic genotyping errors while sensitively detecting genuine recombination transitions.

```bash
python step1.indv.提取个体标记.py -i DB.male.markers.map -g $i -o ./01.indv/${i}.indv.txt
python step2.indv.确定重组位置.py -i 01.indv/${i}.indv.txt -w 2000 -o 02.results/${i}.cross.txt
```
```python
#step1.indv.提取个体标记.py
import argparse
import sys

def extract_columns(input_file, gene_marker, output_file):
    """
    从制表符分隔的表中提取指定列
    
    参数:
    input_file: 输入文件路径
    gene_marker: 基因标记名称
    output_file: 输出文件路径
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"错误: 文件 '{input_file}' 为空")
            return
        
        # 解析标题行
        header = lines[0].strip().split('\t')
        
        # 检查必需的列是否存在
        required_columns = ['#CHROM', 'POS', gene_marker]
        missing_columns = []
        
        for col in required_columns:
            if col not in header:
                missing_columns.append(col)
        
        if missing_columns:
            print(f"错误: 表中缺少以下列: {', '.join(missing_columns)}")
            print(f"可用列: {', '.join(header)}")
            return
        
        # 获取列索引
        chrom_idx = header.index('#CHROM')
        pos_idx = header.index('POS')
        gene_idx = header.index(gene_marker)
        
        # 写入输出文件
        with open(output_file, 'w') as f:
            # 写入标题行
            f.write(f"#CHROM\tPOS\t{gene_marker}\n")
            
            # 写入数据行
            for i, line in enumerate(lines[1:], 1):
                try:
                    fields = line.strip().split('\t')
                    if len(fields) <= max(chrom_idx, pos_idx, gene_idx):
                        print(f"警告: 第{i}行列数不足，跳过此行")
                        continue
                    
                    chrom = fields[chrom_idx]
                    pos = fields[pos_idx]
                    gene_value = fields[gene_idx]
                    
                    f.write(f"{chrom}\t{pos}\t{gene_value}\n")
                    
                except Exception as e:
                    print(f"警告: 处理第{i}行时出错: {e}，跳过此行")
        
        print(f"成功提取数据到: {output_file}")
        print(f"提取的列: #CHROM, POS, {gene_marker}")
        
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{input_file}'")
    except Exception as e:
        print(f"错误: 处理文件时出错: {e}")

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(
        description='从制表符分隔的表中提取指定基因标记的数据',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python work.py -i data.txt -g PD2-101 -o PD2-101.txt
  python work.py -i data.txt -g PD2-10 -o PD2-10.txt
        """
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='输入文件路径（制表符分隔）')
    parser.add_argument('-g', '--gene', required=True, 
                       help='基因标记名称（如: PD2-101, PD2-10）')
    parser.add_argument('-o', '--output', required=True, 
                       help='输出文件路径')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用提取函数
    extract_columns(args.input, args.gene, args.output)

if __name__ == "__main__":
    main()
```

```python
#step2.indv.确定重组位置.py
#!/usr/bin/env python3
"""
单倍型转换点检测脚本
按染色体检测单倍型转换点，支持命令行参数
使用方法: python work.py -i input.txt -o output.txt
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
from collections import defaultdict

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='检测单倍型转换点（按染色体分析）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python work.py -i input.txt -o results.txt
  python work.py -i input.txt -o results.txt -w 200 -t1 0.5 -t2 0.5
  python work.py -i input.txt -o results.txt --window 100 --min_windows 2
  
输入文件格式（制表符分隔，带表头）:
  chr     pos     sample_name
  chr1    1000    DBhap1
  chr1    2000    DBhap1
  chr1    3000    DBhap2
        """
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='输入文件路径（必须，三列：chr, pos, hap）')
    parser.add_argument('-o', '--output', required=True,
                       help='输出文件路径（必须）')
    parser.add_argument('-w', '--window', type=int, default=2000,
                       help='滑动窗口大小（默认: 2000）')
    parser.add_argument('-t1', '--threshold_hap1', type=float, default=0.5,
                       help='hap1主导阈值，hap2比例小于此值判为hap1区（默认: 0.5）')
    parser.add_argument('-t2', '--threshold_hap2', type=float, default=0.5,
                       help='hap2主导阈值，hap2比例大于此值判为hap2区（默认: 0.5）')
    parser.add_argument('-m', '--min_windows', type=int, default=3,
                       help='最小连续窗口数（默认: 3）')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='显示详细处理信息')
    
    return parser.parse_args()

def detect_hap_transitions(df, sample_name, window_size=2000, 
                          hap2_threshold=0.5, hap1_threshold=0.5, 
                          min_consensus_windows=3, verbose=False):
    """
    按染色体检测单倍型转换点
    
    参数:
    ----------
    df : pandas.DataFrame
        包含三列的数据框: chr, pos, hap
    sample_name : str
        样本名
    window_size : int, 默认 2000
        滑动窗口大小
    hap2_threshold : float, 默认 0.5
        hap2主导阈值
    hap1_threshold : float, 默认 0.5
        hap1主导阈值
    min_consensus_windows : int, 默认 3
        最小连续窗口数
    verbose : bool, 默认 False
        是否显示详细信息
    
    返回:
    ----------
    list : 转换点信息列表
    """
    
    # 获取列名
    chr_col = df.columns[0]
    pos_col = df.columns[1]
    hap_col = df.columns[2]
    
    all_transitions = []
    
    for chr_name, chr_df in df.groupby(chr_col):
        if verbose:
            print(f"  分析染色体: {chr_name}，数据点数: {len(chr_df)}", end="")
        
        # 对每个染色体单独分析
        indices = chr_df[pos_col].values
        haps = chr_df[hap_col].values
        
        # 转换为数值（支持DBhap1和DBhap2格式）
        hap_numeric = []
        valid_count = 0
        
        for h in haps:
            if isinstance(h, str):
                h_lower = h.lower()
                if 'hap1' in h_lower or h == '1' or h == '0':
                    hap_numeric.append(0)
                    valid_count += 1
                elif 'hap2' in h_lower or h == '2':
                    hap_numeric.append(1)
                    valid_count += 1
                else:
                    hap_numeric.append(np.nan)
            else:
                if pd.isna(h):
                    hap_numeric.append(np.nan)
                elif h == 0 or h == 1:
                    hap_numeric.append(0)
                    valid_count += 1
                elif h == 2:
                    hap_numeric.append(1)
                    valid_count += 1
                else:
                    hap_numeric.append(np.nan)
        
        hap_numeric = np.array(hap_numeric)
        
        # 检查是否有足够的数据
        if valid_count < window_size:
            if verbose:
                print(f" → 跳过（有效点{valid_count} < 窗口{window_size}）")
            continue
        
        # 滑动窗口计算hap2比例
        half_window = window_size // 2
        n = len(hap_numeric)
        hap2_ratios = np.full(n, np.nan)
        
        for i in range(n):
            start = max(0, i - half_window)
            end = min(n, i + half_window + 1)
            window_vals = hap_numeric[start:end]
            non_nan_vals = window_vals[~np.isnan(window_vals)]
            
            if len(non_nan_vals) > 0:
                hap2_ratios[i] = np.sum(non_nan_vals) / len(non_nan_vals)
        
        # 根据阈值判定状态
        states = np.full(n, -1, dtype=int)  # -1:混合, 0:hap1, 1:hap2
        
        for i in range(n):
            if not np.isnan(hap2_ratios[i]):
                if hap2_ratios[i] > hap2_threshold:
                    states[i] = 1
                elif hap2_ratios[i] < hap1_threshold:
                    states[i] = 0
        
        # 识别连续状态区域
        regions = []
        
        # 找到第一个非混合状态
        current_state = None
        region_start = 0
        for i in range(n):
            if states[i] != -1:
                current_state = states[i]
                region_start = i
                break
        
        if current_state is None:
            if verbose:
                print(" → 跳过（无明确状态）")
            continue
        
        for i in range(region_start + 1, n):
            if states[i] != current_state:
                if current_state != -1:
                    region_length = i - region_start
                    # 计算区域内有效点的比例
                    valid_in_region = np.sum(~np.isnan(hap_numeric[region_start:i]))
                    if valid_in_region >= min_consensus_windows:
                        regions.append({
                            'start_idx': region_start,
                            'end_idx': i - 1,
                            'state': current_state,
                            'length': region_length
                        })
                current_state = states[i]
                region_start = i
        
        # 最后一个区域
        if current_state != -1 and region_start < n:
            region_length = n - region_start
            valid_in_region = np.sum(~np.isnan(hap_numeric[region_start:n]))
            if valid_in_region >= min_consensus_windows:
                regions.append({
                    'start_idx': region_start,
                    'end_idx': n - 1,
                    'state': current_state,
                    'length': region_length
                })
        
        # 检测转换点
        if len(regions) < 2:
            if verbose:
                print(" → 无转换点")
            continue
        
        chr_transitions = []
        
        for i in range(1, len(regions)):
            prev_region = regions[i-1]
            curr_region = regions[i]
            
            if prev_region['state'] == 0 and curr_region['state'] == 1:
                transition_type = 'hap1->hap2'
            elif prev_region['state'] == 1 and curr_region['state'] == 0:
                transition_type = 'hap2->hap1'
            else:
                continue
            
            # 转换位置
            trans_start_idx = prev_region['end_idx']
            trans_end_idx = curr_region['start_idx']
            
            # 获取对应的原始位置
            trans_start_pos = int(indices[trans_start_idx])
            trans_end_pos = int(indices[trans_end_idx])
            
            # 获取区域的实际位置范围
            hap1_start = int(indices[prev_region['start_idx']]) if prev_region['state'] == 0 else int(indices[curr_region['start_idx']])
            hap1_end = int(indices[prev_region['end_idx']]) if prev_region['state'] == 0 else int(indices[curr_region['end_idx']])
            hap2_start = int(indices[prev_region['start_idx']]) if prev_region['state'] == 1 else int(indices[curr_region['start_idx']])
            hap2_end = int(indices[prev_region['end_idx']]) if prev_region['state'] == 1 else int(indices[curr_region['end_idx']])
            
            transition_info = {
                'sample': sample_name,
                'chr': str(chr_name),
                'type': transition_type,
                'transition_start': trans_start_pos,
                'transition_end': trans_end_pos,
                'hap1_start': hap1_start,
                'hap1_end': hap1_end,
                'hap2_start': hap2_start,
                'hap2_end': hap2_end,
                'hap1_length': prev_region['length'] if prev_region['state'] == 0 else curr_region['length'],
                'hap2_length': prev_region['length'] if prev_region['state'] == 1 else curr_region['length']
            }
            
            chr_transitions.append(transition_info)
        
        all_transitions.extend(chr_transitions)
        
        if verbose:
            print(f" → 找到 {len(chr_transitions)} 个转换点")
    
    return all_transitions

def save_results(transitions, output_file, verbose=False):
    """保存结果到文件"""
    if not transitions:
        if verbose:
            print("警告: 未检测到任何转换点，创建空输出文件")
        # 创建空文件
        with open(output_file, 'w') as f:
            f.write("sample\tchr\ttype\ttransition_start\ttransition_end\thap1_start\thap1_end\thap2_start\thap2_end\thap1_length\thap2_length\n")
        return
    
    with open(output_file, 'w') as f:
        # 写文件头
        f.write("sample\tchr\ttype\ttransition_start\ttransition_end\t")
        f.write("hap1_start\thap1_end\thap2_start\thap2_end\t")
        f.write("hap1_length\thap2_length\n")
        
        # 写每个转换
        for trans in transitions:
            f.write(f"{trans['sample']}\t{trans['chr']}\t{trans['type']}\t")
            f.write(f"{trans['transition_start']}\t{trans['transition_end']}\t")
            f.write(f"{trans['hap1_start']}\t{trans['hap1_end']}\t")
            f.write(f"{trans['hap2_start']}\t{trans['hap2_end']}\t")
            f.write(f"{trans['hap1_length']}\t{trans['hap2_length']}\n")
    
    if verbose:
        print(f"结果已保存到: {output_file}")

def main():
    """主函数"""
    # 解析命令行参数
    args = parse_arguments()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在")
        sys.exit(1)
    
    # 检查参数有效性
    if args.window <= 0:
        print("错误: 窗口大小必须为正整数")
        sys.exit(1)
    
    if args.threshold_hap1 < 0 or args.threshold_hap1 > 1 or \
       args.threshold_hap2 < 0 or args.threshold_hap2 > 1:
        print("错误: 阈值必须在0到1之间")
        sys.exit(1)
    
    if args.threshold_hap1 >= args.threshold_hap2:
        print("警告: hap1阈值应小于hap2阈值，否则可能无明确区域")
    
    if args.min_windows < 1:
        print("错误: 最小窗口数必须为正整数")
        sys.exit(1)
    
    # 显示参数信息
    print("=" * 70)
    print("单倍型转换点检测")
    print("=" * 70)
    print(f"输入文件: {args.input}")
    print(f"输出文件: {args.output}")
    print(f"参数设置:")
    print(f"  - 窗口大小: {args.window}")
    print(f"  - hap1阈值: <{args.threshold_hap1}")
    print(f"  - hap2阈值: >{args.threshold_hap2}")
    print(f"  - 最小连续窗口: {args.min_windows}")
    print("-" * 70)
    
    try:
        # 读取数据
        if args.verbose:
            print("正在读取数据...")
        
        df = pd.read_csv(args.input, sep='\t')
        
        if len(df.columns) < 3:
            print(f"错误: 输入文件需要至少3列，当前有{len(df.columns)}列")
            print("列名:", list(df.columns))
            sys.exit(1)
        
        # 获取样本名（第三列表头）
        sample_name = df.columns[2]
        
        if args.verbose:
            print(f"样本: {sample_name}")
            print(f"总数据点数: {len(df)}")
            print(f"染色体数: {df.iloc[:, 0].nunique()}")
            print("开始检测转换点...")
        
        # 检测转换点
        transitions = detect_hap_transitions(
            df=df,
            sample_name=sample_name,
            window_size=args.window,
            hap2_threshold=args.threshold_hap2,
            hap1_threshold=args.threshold_hap1,
            min_consensus_windows=args.min_windows,
            verbose=args.verbose
        )
        
        # 保存结果
        save_results(transitions, args.output, args.verbose)
        
        # 汇总统计
        print("\n" + "=" * 70)
        print("检测完成!")
        print("=" * 70)
        
        if transitions:
            # 按染色体统计
            chr_stats = defaultdict(int)
            type_stats = defaultdict(int)
            
            for trans in transitions:
                chr_stats[trans['chr']] += 1
                type_stats[trans['type']] += 1
            
            print(f"总共检测到 {len(transitions)} 个转换点")
            print(f"转换类型统计:")
            for trans_type, count in type_stats.items():
                print(f"  {trans_type}: {count} 个")
            
            print(f"\n染色体统计:")
            for chr_name, count in sorted(chr_stats.items()):
                print(f"  {chr_name}: {count} 个转换点")
            
            # 显示前几个转换点
            if args.verbose and len(transitions) > 0:
                print(f"\n前5个转换点:")
                for i, trans in enumerate(transitions[:5], 1):
                    print(f"{i}. {trans['sample']} {trans['chr']} {trans['type']}:")
                    print(f"   位置: {trans['transition_start']:,} - {trans['transition_end']:,}")
                    print(f"   hap1: {trans['hap1_start']:,}-{trans['hap1_end']:,} (长度: {trans['hap1_length']})")
                    print(f"   hap2: {trans['hap2_start']:,}-{trans['hap2_end']:,} (长度: {trans['hap2_length']})")
            
            print(f"\n详细结果已保存到: {args.output}")
        else:
            print("未检测到任何转换点")
            print("建议: 尝试调整参数（减小窗口大小或调整阈值）")
        
        print("=" * 70)
        
    except FileNotFoundError as e:
        print(f"文件错误: {e}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("错误: 输入文件为空")
        sys.exit(1)
    except pd.errors.ParserError as e:
        print(f"解析错误: {e}")
        print("请确保输入文件是制表符分隔的文本文件")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行错误: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
```
```bash
cat PD2-290.cross.txt 
sample	chr	type	transition_start	transition_end	hap1_start	hap1_end	hap2_start	hap2_end	hap1_length	hap2_length
PD2-290	chr1	hap1->hap2	14938215	14938365	46150	14938215	14938365	32200565	72469	57223
PD2-290	chr2	hap1->hap2	4633066	4633350	214565	4633066	4633350	33358895	21493	133812
PD2-290	chr4	hap2->hap1	372038	372294	372294	30561974	64985	372038	123220	2358
PD2-290	chr4	hap1->hap2	30561974	30562290	372294	30561974	30562290	32858183	123220	13721
PD2-290	chr5	hap2->hap1	39136606	39136856	39136856	48916282	307533	39136606	47569	176941
PD2-290	chr6	hap1->hap2	20448408	20448456	11536	20448408	20448456	27080315	81672	35289
PD2-290	chr7	hap1->hap2	89616	89829	84502	89616	89829	17903574	163	69033
PD2-290	chr7	hap2->hap1	17903574	17904276	17904276	31654137	89829	17903574	56961	69033
PD2-290	chr8	hap1->hap2	4401684	4401830	74300	4401684	4401830	21241351	26917	59069
PD2-290	chr8	hap2->hap1	21241351	21241390	21241390	30726327	4401830	21241351	46139	59069
PD2-290	chr9	hap1->hap2	8377343	8377570	9255	8377343	8377570	30192576	51648	102365
```
---

## Step 5 — Population Recombination Rate Analysis

### 5.1 Variant Phasing (BEAGLE v5.5)

Accessions from the DB-related (apomictic) and PN-related (sexual) subgroups were phased independently using BEAGLE v5.5 under default parameters. Phased haplotypes are required as input for both SMC++ demographic inference and Pyrho recombination rate estimation.

```bash
java -Xmx500G -jar beagle.22Jul22.46e.jar gt=${SUBGROUP}.vcf.gz out=${SUBGROUP}.phased nthreads=${THREADS}
```

---

### 5.2 Demographic History Inference (SMC++ v1.15.5)

SMC++ was used to infer the effective population size history of each subgroup independently, using a mutation rate of μ = 2.2 × 10⁻⁸ per base per generation consistent with prior *Citrus* studies. The resulting demographic trajectories serve as priors to improve the accuracy of downstream recombination rate estimation with Pyrho.

```bash
python3 -i BEAGLE.vcf -l sample_prefix -o output.file -f genome.fa.fai

```
```python
#smc-analysis.py
import vcf
import os
import argparse

parser = argparse.ArgumentParser(description = 'For citrus', add_help = False, usage = '\npython3 -i [input file from step1] -l [sample list] -o [output.file]')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i', '--input', metavar = '[input.beagle.vcf]', help = 'input.beagle.vcf', required = True)
required.add_argument('-l', '--fam', metavar = '[fam name]', help = 'fam_name', required = True)
required.add_argument('-f', '--fai', metavar = '[refence.fai]', help = 'refence.fai', required = True)
required.add_argument('-m', '--mu', metavar = '[mutation ratio]', help = 'mutation ratio', default='2.2e-8',type=str)
optional.add_argument('-h', '--help', action = 'help', help = 'help')

args = parser.parse_args()

#os.system()
family=args.fam
sample_list=[]
mu=args.mu
vcf_reader = vcf.Reader(open(args.input, 'r'))
for record in vcf_reader:
    for sample in record.samples:
        #print(sample.sample)
        sample_list.append(sample.sample)
    break

dict={}
with open(args.fai,'r') as fai_file:
    for i in fai_file:
        i=i.replace("\n","").split('\t')
        chr=i[0]
        changdu=int(i[1])
        dict[chr]=changdu

cmd3='mkdir $PWD/{}_out'.format(family)
#cmd2='docker pull docker-0.unsee.tech/terhorst/smcpp'
#print(cmd2)
#os.system(cmd2)
#print(cmd3)
os.system(cmd3)
for chr in dict.keys():
    out_tmp1=args.input+chr+'.vcf'
    out_file1=open(out_tmp1,'w')
    header_chr="##contig=<ID="+chr+",length="+str(dict[chr])+">"+'\n'
    header0="##fileformat=VCFv4.2"+'\n'
    out_file1.write(header0)
    out_file1.write(header_chr)
    with open(args.input,'r') as vcf:
        for line in vcf:
            if "#CHROM" in line or line.split('\t')[0]==chr:
                out_file1.write(line)
    cmd0='bgzip -c {} > {}.gz'.format(out_tmp1,out_tmp1)
    #print(cmd0)
    os.system(cmd0)
    cmd1='tabix -p vcf {}.gz'.format(out_tmp1)
    #print(cmd1)
    os.system(cmd1)
    
    s=''
    for one in sample_list:
        s+=one+','
    s=s[:-1]
    cmd4='docker run -v $PWD:/work -w /work terhorst/smcpp:latest vcf2smc /work/{}.gz {}_out/{}.smc.gz {} {}:{}'.format(out_tmp1,family,out_tmp1,chr,family,s)
    #print(cmd4)
    os.system(cmd4)
    
cmd5='mkdir $PWD/{}_analysis'.format(family)
#print(cmd5)
os.system(cmd5)
##for single chr
s1=''
s2=''
for chr in dict.keys():
    s1+='/work/'+family+'_out/'+args.input+chr+'.vcf'+'.smc.gz'+' '
    s2+='/work/{}_analysis/{}/model.final.json'.format(family,chr)+' '
    out_tmp1=args.input+chr+'.vcf'
    cmd6='docker run -v $PWD:/work -w /work terhorst/smcpp:latest estimate -o /work/{}_analysis/{} {} /work/{}_out/{}.smc.gz'.format(family,chr,mu,family,out_tmp1)
    #print(cmd6)
    os.system(cmd6)
##for combine
cmd7='mkdir $PWD/{}_combine'.format(family)
#print(cmd7)
os.system(cmd7)
s1=s1[:-1]
s2=s2[:-1]
cmd8='docker run -v $PWD:/work -w /work terhorst/smcpp:latest estimate -o /work/{}_combine/ {} {}'.format(family,mu,s1)
#print(cmd8)
os.system(cmd8)
cmd9='docker run -v $PWD:/work -w /work terhorst/smcpp:latest plot /work/{}.pdf /work/{}_combine/model.final.json {}'.format(family,family,s2)
#print(cmd9)
os.system(cmd9)
cmd10='rm {}*.vcf*'.format(family)
#os.system(cmd10)
```


---

### 5.3 Population Recombination Rate Estimation (Pyrho)

Fine-scale population recombination rates were estimated for the DB-related and PN-related subgroups using Pyrho, with approximate lookup tables constructed from demographic history. For each subgroup, the downsample size was set to n = 28 haplotypes and the lookup table size to N = 50, balancing computational efficiency and estimation accuracy. Recombination rates were estimated independently for each of the nine chromosomes.

```bash
singularity run /home/wangnan/Software/pyrho_20251120.sif pyrho make_table -n 34 -N 42 --mu 2.2e-8 --logfile /home/wangnan/8.kumquat/10.natural.population/03.pyrho/logs_PN_34_42.txt --outfile PN_n_34_N_42_2026_01_12_lookuptable.hdf --approx --numthreads 10 --smcpp_file PN-plot_chr2.csv --decimate_rel_tol 0.1
singularity run /home/wangnan/Software/pyrho_20251120.sif pyrho make_table -n 40 -N 50 --mu 2.2e-8 --logfile /home/wangnan/8.kumquat/10.natural.population/03.pyrho/logs_APO_40_50.txt --outfile APO_n_40_N_50_2026_01_12_lookuptable.hdf --approx --numthreads 10 --smcpp_file APO-plot_chr6.csv --decimate_rel_tol 0.1

for i in {1..9};do singularity run /home/wangnan/Software/pyrho_20251120.sif pyrho optimize --vcffile /home/wangnan/8.kumquat/10.natural.population/03.pyrho/PN.chr${i}.vcf --windowsize 50 --blockpenalty 50 --tablefile PN_n_34_N_42_2026_01_12_lookuptable.hdf --ploidy 1 --outfile PN_34_chr${i}.Recom.txt --numthreads 10;done
for i in {1..9};do singularity run /home/wangnan/Software/pyrho_20251120.sif pyrho optimize --vcffile /home/wangnan/8.kumquat/10.natural.population/03.pyrho/APO.chr${i}.vcf --windowsize 50 --blockpenalty 50 --tablefile APO_n_40_N_50_2026_01_12_lookuptable.hdf --ploidy 1 --outfile APO_40_chr${i}.Recom.txt --numthreads 10;done

```

---

## Step 6 — Haplotype Shuffling Map and Segregation Distortion Analysis

### 6.1 Genomic Bin Construction

The genome was partitioned into non-overlapping bins using recombination breakpoints as boundaries, such that each bin represents a segment free of internal crossover events. Within each bin, the haplotype combination of each F1 offspring was classified into one of four states according to the parental haplotype combination (PN_hap1/DB_hap1, PN_hap1/DB_hap2, PN_hap2/DB_hap1, PN_hap2/DB_hap2), forming the haplotype shuffling map.

```bash
python 提取构建洗牌关系的脚本.py [-h] -i/--input, -c/--cross, -o/--output

OUTPUT
cat 03.final_maps/PD2-325.haps.STAT 
Sample	Chr	Start	End	Hap1_Ratio	Hap2_Ratio
./01.indv/PD2-325.indv.txt	chr1	0	3198232	0.996581	0.001681
./01.indv/PD2-325.indv.txt	chr1	3198232	14454951	0.002280	0.993011
./01.indv/PD2-325.indv.txt	chr1	14454951	32266151	0.990246	0.002746
./01.indv/PD2-325.indv.txt	chr2	0	9125915	0.991059	0.002845
./01.indv/PD2-325.indv.txt	chr2	9125915	33368013	0.002504	0.995046
./01.indv/PD2-325.indv.txt	chr3	0	41283659	0.989891	0.004513
./01.indv/PD2-325.indv.txt	chr4	0	32876879	0.994953	0.002656
./01.indv/PD2-325.indv.txt	chr5	0	10329925	0.002166	0.995846
./01.indv/PD2-325.indv.txt	chr5	10329925	43408632	0.992700	0.003942
./01.indv/PD2-325.indv.txt	chr5	43408632	48982865	0.002591	0.994779
./01.indv/PD2-325.indv.txt	chr6	0	22226616	0.002803	0.993069
./01.indv/PD2-325.indv.txt	chr6	22226616	27141800	0.994943	0.001686
./01.indv/PD2-325.indv.txt	chr7	0	395736	0.986315	0.013685
./01.indv/PD2-325.indv.txt	chr7	395736	31662302	0.002923	0.991360
./01.indv/PD2-325.indv.txt	chr8	0	30741539	0.990978	0.004564
./01.indv/PD2-325.indv.txt	chr9	0	3996260	0.996668	0.001794
./01.indv/PD2-325.indv.txt	chr9	3996260	30315542	0.002707	0.993031

```
```python
#提取构建洗牌关系的脚本.py
#!/usr/bin/env python3
import argparse
import sys

# 染色体长度定义
Chr_len = {
    'chr1': 32266151, 'chr2': 33368013, 'chr3': 41283659,
    'chr4': 32876879, 'chr5': 48982865, 'chr6': 27141800,
    'chr7': 31662302, 'chr8': 30741539, 'chr9': 30315542
}

# 解析命令行参数
parser = argparse.ArgumentParser(description='处理染色体数据')
parser.add_argument('-i', '--input', required=True, help='输入文件路径')
parser.add_argument('-c', '--cross', required=True, help='交叉文件路径')
parser.add_argument('-o', '--output', required=True, help='输出文件路径')
args = parser.parse_args()

indv_file = args.input
cross_file = args.cross
output_file = args.output

def filter_chr_region(file_path, chrom, start, end):
    """筛选指定染色体区域的数据"""
    result_list = []
    hap1_count = 0
    hap2_count = 0
    
    with open(file_path, 'r') as f:
        for line in f:
            if not line.strip():
                continue
                
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            
            chr_col, pos_col, value_col = parts[0], parts[1], parts[2]
            
            if chr_col == chrom:
                position = int(pos_col)
                if start <= position <= end:
                    result_list.append(value_col)
                    if "hap1" in value_col:
                        hap1_count += 1
                    elif "hap2" in value_col:
                        hap2_count += 1
    
    # 计算比例
    if not result_list:
        return 0.0, 0.0
    
    total = len(result_list)
    hap1_ratio = hap1_count / total
    hap2_ratio = hap2_count / total
    
    return hap1_ratio, hap2_ratio

# 打开输出文件
outputFile = open(output_file, 'w')
outputFile.write('Sample\tChr\tStart\tEnd\tHap1_Ratio\tHap2_Ratio\n')

# 处理每个染色体
for i in range(1, 10):
    Chr = 'chr' + str(i)
    
    # 读取交叉文件，获取位置列表
    list_pos = [0]
    with open(cross_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
                
            line_data = line.strip().split('\t')
            if len(line_data) < 5:
                continue
                
            if line_data[1] == Chr:
                pos1 = int(line_data[3])
                pos2 = int(line_data[4])
                mid_point = 0.5 * (pos1 + pos2)
                list_pos.append(int(mid_point))
    
    # 添加染色体末端
    list_pos.append(Chr_len[Chr])
    list_pos = sorted(set(list_pos))  # 去重并排序
    
    # 处理每个区间
    for j in range(len(list_pos) - 1):
        start_pos = list_pos[j]
        end_pos = list_pos[j + 1]
        
        if end_pos > start_pos:  # 跳过无效区间
            hap1_ratio, hap2_ratio = filter_chr_region(indv_file, Chr, start_pos, end_pos)
            
            # 写入结果
            result_line = f"{indv_file}\t{Chr}\t{start_pos}\t{end_pos}\t{hap1_ratio:.6f}\t{hap2_ratio:.6f}"
            outputFile.write(result_line + '\n')

outputFile.close()
print(f"处理完成！结果保存到: {output_file}")
```

---

## Step 7 — Mapping of Causal Genetic Factors

### 7.1 Lethal Interval Identification

The genome-wide distortion index (differential allelic frequency between surviving and non-surviving offspring) was scanned across all bins to identify the region with maximum distortion as the candidate lethal interval. Combined with haplotype transmission patterns, DB_hap2 was confirmed to carry the causal lethal variation. The minimal lethal interval was then precisely defined based on the surrounding breakpoint coordinates.

```bash
cd /home/wangnan/8.kumquat/13.SIFT_based_deleterious/1.PlantCAD2
python arg.循环版step2.grep.specific.variants_副本.py -i $i -o RESULTS_${i}.TXT

```
```python
#arg.循环版step2.grep.specific.variants_副本.py 
import os
import argparse
#解析命令行参数
parser = argparse.ArgumentParser(description='查找变异对应的CSV文件')
parser.add_argument('-i', '--input', required=True, help='输入文件路径')
parser.add_argument('-o', '--output', required=True, help='输出文件路径')

args = parser.parse_args()
def csv_find(chrom,position,index_file="Z.duiying.txt"):
    out_file=''
    with open(index_file,'r') as f0:
        for line0 in f0:
            line0=line0.replace("\n","").split('\t')
            #print(line0)
            chr_index="chr"+line0[0]
            position_start=int(line0[1])
            position_end=int(line0[2])
            csv_file=line0[3]
            if chr_index==chrom and position_start<=position and position<=position_end:
                out_file=csv_file
    return out_file
#csv_find('chr1',8852844)

def from_csv_to_dict(chrom,csv_file):
    dict={}
    with open(csv_file,'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            line=line.replace("\n","").split(',')
            #print(line)
            chr=chrom
            position=int(line[0])
            ref=line[1]
            alt=line[2]
            dict[chr+'_'+str(position)+'_'+ref+'_'+alt]=",".join(line[-4:])
    return dict

outfile=open(args.output,'w')
with open(args.input,'r') as f2:
    for line in f2:
        line=line.replace("\n","").split('\t')
        chr=line[0]
        start=line[1]
        end=line[2]
        name=line[0]+'_'+line[1]+'_'+line[2]
        d0=float(line[4])
        R2=line[5]
        if d0>0:#测试hap2的有害变异
            test_vcf="vcftools --vcf merged_DBhap2.vcf --recode --stdout --chr "+chr+" --from-bp "+start+" --to-bp "+end+" > "+name+"_DBhap2.vcf"
            #print(test_vcf)
            os.system(test_vcf)
            #with open(name+"_DBhap2.vcf",'r') as f3:
            with open(name+"_DBhap2.vcf",'r') as f3:
                lines3=f3.readlines()
                for line3 in lines3:
                    if not line3.startswith("#"):
                        line3=line3.replace("\n","").split('\t')
                        #print(line3)
                        chr=line3[0]
                        position=int(line3[1])
                        ref=line3[3]
                        alt=line3[4]
                        csv_file=csv_find(chr,position)
                        #print(csv_file)
                        csv_dict=from_csv_to_dict(chr,"./database/"+csv_file)
                        if chr+'_'+str(position)+'_'+ref+'_'+alt in csv_dict.keys():
                            #print(csv_dict[chr+'_'+str(position)+'_'+ref+'_'+alt])
                            line3.append(csv_dict[chr+'_'+str(position)+'_'+ref+'_'+alt])
                        else:
                            #print("Notfound")
                            line3.append("Notfound")
                        outfile.write("\t".join(str(i) for i in line3)+"\n")
                        outfile.flush()
            os.system("rm "+name+"_DBhap2.vcf")
        if d0<0:#测试hap1的有害变异
            test_vcf="vcftools --vcf merged_DBhap1.vcf --recode --stdout --chr "+chr+" --from-bp "+start+" --to-bp "+end+" > "+name+"_DBhap1.vcf"
            #print(test_vcf)
            os.system(test_vcf)
            with open(name+"_DBhap1.vcf",'r') as f3:
                lines3=f3.readlines()
                for line3 in lines3:
                    if not line3.startswith("#"):
                        line3=line3.replace("\n","").split('\t')
                        #print(line3)
                        chr=line3[0]
                        position=int(line3[1])
                        ref=line3[3]
                        alt=line3[4]
                        csv_file=csv_find(chr,position)
                        #print(csv_file)
                        csv_dict=from_csv_to_dict(chr,"./database/"+csv_file)
                        if chr+'_'+str(position)+'_'+ref+'_'+alt in csv_dict.keys():
                            #print(csv_dict[chr+'_'+str(position)+'_'+ref+'_'+alt])
                            line3.append(csv_dict[chr+'_'+str(position)+'_'+ref+'_'+alt])
                        else:
                            #print("Notfound")
                            line3.append("Notfound")
                        outfile.write("\t".join(str(i) for i in line3)+"\n")
                        outfile.flush()
            os.system("rm "+name+"_DBhap1.vcf")
```

---

### 7.3 Functional Effect Scoring — SIFT4G

SIFT4G predicts variant deleteriousness by measuring amino acid conservation across homologous sequences: variants at highly conserved positions (SIFT score < 0.05) are classified as deleterious. A species-specific SIFT4G database was first constructed from the *C. hindsii* proteome and a curated plant protein database, then applied to all genic candidate variants in the lethal interval.

```bash

# Annotate variants
java -jar SIFT4G_Annotator.jar -c -i DB_hap2.lethal_genic.vcf.gz -d ${SIFT4G_DB_DIR} -r sift4g_results/
```

---

### 7.4 Functional Effect Scoring — PlantCAD2 (Zero-shot LLM Scoring)

PlantCAD2 is a plant-specific genomic large language model that performs zero-shot variant pathogenicity prediction by leveraging cross-species sequence conservation and genomic context features, without requiring task-specific training. All candidate SNPs in the lethal interval were scored to provide a complementary pathogenicity assessment independent of the alignment-based SIFT approach.

```bash
python plantcad_score.py --vcf DB_hap2.lethal_genic.vcf.gz --genome DB_hap2.final.fa --model plantcad2_pretrained/ --out plantcad2_scores.tsv
```

---
