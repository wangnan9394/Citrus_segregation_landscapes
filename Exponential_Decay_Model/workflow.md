# Exponential Decay Model of Segregation Distortion

## Part I — Single-locus Exponential Decay Model

### Model Rationale

We developed a single-locus exponential decay model to describe the spatial topography of transmission distortion surrounding a causal segregating locus. The core principle is that while a distorting allele (e.g., a gametic lethal) exerts a primary selective pressure at its physical position, its influence on the surrounding genomic landscape is mediated by the **local recombination frequency** — neutral markers farther from the causal site are more likely to be decoupled by recombination.

### Distortion Index Definition

The distortion index *D* at genomic position *x* relative to a distorting locus is defined as the normalized difference between the favored haplotype frequency *f₁(x)* and the unfavored haplotype frequency *f₂(x)*:

$$D_x = \frac{f_1(x) - f_2(x)}{f_1(x) + f_2(x)}$$

At the causal site (*x* = *POS*), the distortion reaches its maximum value *D₀*, representing the raw selective coefficient. As physical distance increases, recombination decouples markers from the causal locus, leading to exponential decay in *D*.

### Decay Equation

$$D_x = D_0 \times e^{-r \times \frac{|x - \text{POS}|}{100}}$$

| Symbol | Meaning | Units |
|--------|---------|-------|
| *D_x* | Distortion index at position *x* | dimensionless [−1, 1] |
| *D₀* | Maximum distortion at the causal site (selective strength) | dimensionless [−1, 1] |
| *r* | Local recombination rate | cM / Mb |
| *x* | Query genomic position | cM |
| POS | Physical position of the causal locus | cM |
| \|x − POS\|| Absolute distance from causal locus | cM |
| 100 | Scaling factor for distance normalization | dimensionless |

The scaling factor `/100` converts the product `r · distance` into a dimensionless exponent, ensuring that decay occurs over biologically realistic distances given typical plant recombination rates (~1–5 cM/Mb).

```python
# single-locus exponential decay and plots
import numpy as np
import matplotlib.pyplot as plt

# 设置绘图风格
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['Arial']

# 1. 定义物理距离范围 (例如 0-200 Mb, 扩大范围以观察缩放后的衰减)
# 加上 /100 后，衰减变慢，因此我们需要更长的距离来观察趋势
physical_distance = np.linspace(0, 200, 500) 

# 2. 定义不同的重组率参数 r (cM/Mb)
r_values = [0.5, 1.0, 2.0, 5.0]

# 3. 定义初始选择强度 D0
D0_values = [0.3, 0.5, 0.7, 0.9]

# ==========================================
# Figure 1: 不同初始强度下的 D 指数衰减 (含 /100 缩放)
# ==========================================
plt.figure(figsize=(14, 10))

for i, D0 in enumerate(D0_values):
    plt.subplot(2, 2, i+1)
    
    for r in r_values:
        # 严格执行带缩放的公式: Dx = D0 * exp(-r * dist / 100)
        Dx = D0 * np.exp(-r * physical_distance / 100)
        
        plt.plot(physical_distance, Dx, 
                 label=f'r = {r} cM/Mb', linewidth=2)
    
    plt.xlabel('Physical Distance from Causal Locus (Mb)')
    plt.ylabel('Distortion Index ($D_x$)')
    plt.title(f'Initial Selection Strength $D_0$ = {D0}')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=9)
    plt.ylim(0, 1.05)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.suptitle("Scaled Single-locus Exponential Decay Model ($D_x = D_0 \cdot e^{-r \cdot \Delta x / 100}$)", fontsize=16)

# ==========================================
# Figure 2: 频率演变与 D 指数的关系 (Detail View)
# ==========================================
plt.figure(figsize=(11, 7))

D0_fixed = 0.8  # 初始 D 值
r_fixed = 1.5   # 重组率
# 
# 使用缩放公式计算 Dx
Dx_profile = D0_fixed * np.exp(-r_fixed * physical_distance / 100)

# 频率推导: f1 = (D + 1) / 2
f1_profile = (Dx_profile + 1) / 2
f2_profile = 1 - f1_profile

plt.plot(physical_distance, f1_profile, label='Favored Haplotype Frequency ($f_1$)', linewidth=2.5, color='royalblue')
plt.plot(physical_distance, f2_profile, label='Unfavored Haplotype Frequency ($f_2$)', linewidth=2.5, color='salmon')
plt.plot(physical_distance, Dx_profile, label='Distortion Index ($D_x$)', linewidth=3, color='black', linestyle='--')

plt.axhline(y=1.0, color='red', linestyle=':', alpha=0.5, label='Theoretical Max $D$ (1.0)')
plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.3, label='Equilibrium (0.5)')

plt.xlabel('Physical Distance (Mb)')
plt.ylabel('Frequency / Distortion Index ($D$)')
plt.title(f'Scaled Decay Profile (Initial $D_0$={D0_fixed}, r={r_fixed}, Scale=1/100)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ==========================================
# 打印转换表
# ==========================================
print("\n" + "="*55)
print("Reference Table: D to Frequencies")
print("="*55)
print(f"{'D Index':<15} | {'Favored f1':<15} | {'Unfavored f2':<15}")
print("-" * 55)
for d_val in [0.9, 0.7, 0.5, 0.3, 0.1, 0.0]:
    f1 = (d_val + 1) / 2
    f2 = 1 - f1
    print(f"{d_val:<15.1f} | {f1:<15.3f} | {f2:<15.3f}")
print("="*55)
```

---

## Part II — Multi (two and three) -locus Superposition and Saturation Modeling

### Model Extension

Real segregation distortion landscapes often exhibit complex multi-peak patterns driven by multiple linked causal loci. We extended the single-locus framework into a **multi-locus superposition model** where the observed distortion at any coordinate *x* emerges from the cumulative interference of *k* discrete distorting loci.

### Per-locus Potential

For each locus *i*, located at position POS_i with selective strength *D_{0,i}* and local recombination rate *r_i*:

$$D_i = D_{0,i} \times e^{-r_i \times \frac{|x - \text{POS}_i|}{100}}$$

*D_i* represents the raw potential shift contributed by an individual locus, before accounting for biological constraints.

### Saturation Constraint via tanh Transformation

A simple linear summation Σ D_i could theoretically yield predicted values exceeding the biological boundary [−1, 1] of the distortion index (since haplotype frequencies must sum to unity). We therefore apply a hyperbolic tangent (**tanh**) transformation to enforce saturation:

$$D_{\text{exp}}(x) = \tanh\left(\sum_{i=1}^{k} D_i(x)\right)$$

This formulation captures two key biological behaviors:

1. **Synergistic saturation**: Proximal loci with the same sign of *D₀* approach a plateau as their combined effect reaches the terminal limit.
2. **Antagonistic cancellation**: Proximal loci with opposing signs of *D₀* partially cancel each other out, potentially masking moderate-effect alleles.

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 设置绘图风格
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['Arial']

def multi_locus_superposition_model(genome_length=150, 
                                   locus_positions=[45, 55], 
                                   D0_values=[0.40, 0.35], 
                                   r_values=[1.0, 1.0], 
                                   n_points=1000):
    """
    Dx = tanh( Σ (D0_i * exp(-r_i * |x - POS_i| / 100)) )
    """
    positions = np.linspace(0, genome_length, n_points)
    aggregate_potential = np.zeros(len(positions))
    individual_potentials = []

    # 1. 计算每个位点的线性潜力 (Di) - 包含 100 缩放
    for pos, d0, r in zip(locus_positions, D0_values, r_values):
        dist = np.abs(positions - pos)
        di = d0 * np.exp(-r * dist / 100)
        individual_potentials.append(di)
        aggregate_potential += di

    # 2. 应用 tanh 饱和转换
    D_obs = np.tanh(aggregate_potential)
    
    results = pd.DataFrame({
        'position_cM': positions, 
        'D_obs': D_obs, 
        'Aggregate_Potential': aggregate_potential
    })
    for i, pot in enumerate(individual_potentials):
        results[f'Locus_{i+1}_Potential'] = pot
        
    return results

def plot_combined_results(results, locus_positions, D0_values, r_values, title):
    """在一个图中展示独立潜力、总潜力和饱和后的观测值"""
    plt.figure(figsize=(12, 7))
    
    # 1. 绘制各个位点的独立潜力 (Individual Potentials)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']
    for i in range(len(locus_positions)):
        plt.plot(results['position_cM'], results[f'Locus_{i+1}_Potential'], 
                 '--', lw=1.5, color=colors[i % len(colors)], 
                 label=f'Locus {i+1} Potential ($D_0$={D0_values[i]})')
        plt.axvline(x=locus_positions[i], color=colors[i % len(colors)], 
                    linestyle=':', alpha=0.4)

    # 2. 绘制线性叠加后的总潜力 (Unsaturated Aggregate)
    plt.plot(results['position_cM'], results['Aggregate_Potential'], 
             'k:', lw=1.5, alpha=0.6, label='Theoretical Sum ($\sum D_i$)')
    
    # 3. 绘制最终饱和观测值 (Observed Distortion)
    plt.plot(results['position_cM'], results['D_obs'], 
             color='darkred', lw=3, label='Observed Distortion ($D_{obs} = tanh(\sum D_i)$)')
    plt.fill_between(results['position_cM'], 0, results['D_obs'], color='red', alpha=0.1)
    
    # 绘制生物学边界和平衡轴
    plt.axhline(y=1, color='blue', linestyle='--', alpha=0.2, label='Limit $\pm 1$')
    plt.axhline(y=-1, color='blue', linestyle='--', alpha=0.2)
    plt.axhline(y=0, color='black', lw=0.8)
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Genomic Position (cM / Mb)')
    plt.ylabel('Distortion Index ($D_x$)')
    plt.legend(loc='upper right', fontsize='small', frameon=True)
    plt.grid(True, alpha=0.2)
    plt.ylim(min(-1.1, results['Aggregate_Potential'].min()-0.1), 
             max(1.1, results['Aggregate_Potential'].max()+0.1))
    
    plt.tight_layout()
    plt.show()

def parameter_scan_english():
    """执行 3 个指定的模拟案例"""
    test_cases = [
        # (Positions, D0s, rs, Description)
        ([45, 55], [0.40, 0.35], [1.0, 1.0], "Case 1: Strong Positive Overlap"),
        ([25, 75, 120], [0.3, -0.25, 0.15], [1.0, 1.0, 1.0], "Case 2: Three Loci Mixed (Multi-peak)"),
        ([30, 45], [0.45, -0.3], [1.0, 1.0], "Case 3: Close Opposite Effects (Interference)")
    ]
    
    for pos, d0s, rs, desc in test_cases:
        print(f"\n--- Running: {desc} ---")
        results = multi_locus_superposition_model(locus_positions=pos, 
                                                 D0_values=d0s, 
                                                 r_values=rs)
        plot_combined_results(results, pos, d0s, rs, desc)

if __name__ == "__main__":
    parameter_scan_english()
```

---

## Part III — Genomic Deconvolution and Empirical Model Fitting
inpput with genetic distance and distortion index

```python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import argparse
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

# ==========================================
# 1. 核心数学模型 (对齐论文，含 100 缩放)
# ==========================================
def single_locus_model(X, D0, r, pos):
    return D0 * np.exp(-r * np.abs(X - pos) / 100)

def additive_model_k_loci(X, *params):
    k = len(params) // 3
    aggregate_potential = np.zeros(len(X), dtype=np.float64)
    for i in range(k):
        D0, r, pos = params[3*i : 3*i+3]
        aggregate_potential += single_locus_model(X, D0, r, pos)
    return np.tanh(aggregate_potential)

# ==========================================
# 2. 核心算法：递归步进搜索 (含 100 缩放)
# ==========================================
def find_loci_stepwise(X, Y, max_k=10, dist_threshold=8):
    current_Y = Y.copy()
    found_params = []
    x_min, x_max = X.min(), X.max()
    
    for k in range(max_k):
        peaks_p, _ = find_peaks(current_Y, height=0.015) 
        peaks_n, _ = find_peaks(-current_Y, height=0.015)
        all_candidates = np.concatenate([peaks_p, peaks_n])
        
        if len(all_candidates) == 0: break
            
        best_idx = all_candidates[np.argmax(np.abs(current_Y[all_candidates]))]
        pos_guess = X[best_idx]
        d0_guess = current_Y[best_idx]
        
        if any(abs(pos_guess - p[2]) < dist_threshold for p in found_params):
            mask = (X > pos_guess - dist_threshold) & (X < pos_guess + dist_threshold)
            current_Y[mask] = 0
            continue
            
        try:
            popt, _ = curve_fit(
                single_locus_model, X, current_Y, 
                p0=[d0_guess, 1.0, pos_guess],
                bounds=([-1.5, 0.01, x_min], [1.5, 15.0, x_max])
            )
            found_params.append(popt)
            current_Y -= single_locus_model(X, *popt)
        except:
            break
            
    return np.array(found_params).flatten()

# ==========================================
# 3. 批量处理逻辑
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Multi-chromosome SD Deconvolution with /100 Scaling")
    parser.add_argument("input", nargs="+", help="Input .txt or .csv files")
    parser.add_argument("-o", "--outputdir", default="output_analysis", help="Output directory")
    args = parser.parse_args()

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    all_results = []
    summary_data = []

    for file_path in args.input:
        try:
            filename = os.path.basename(file_path)
            df = pd.read_csv(file_path, sep=None, engine='python')
            X, Y = df.iloc[:, 0].values, df.iloc[:, 1].values
            
            init_p = find_loci_stepwise(X, Y)
            if len(init_p) > 0:
                k_loci = len(init_p) // 3
                lower = ([-1.5, 0.01, X.min()] * k_loci)
                upper = ([1.5, 15.0, X.max()] * k_loci)
                try:
                    final_popt, _ = curve_fit(additive_model_k_loci, X, Y, p0=init_p, bounds=(lower, upper))
                except:
                    final_popt = init_p
            else:
                final_popt = []

            r2 = 0
            if len(final_popt) > 0:
                y_pred = additive_model_k_loci(X, *final_popt)
                r2 = 1 - np.sum((Y - y_pred)**2) / np.sum((Y - np.mean(Y))**2)

            all_results.append((X, Y, final_popt, filename, r2))

            for j in range(len(final_popt)//3):
                summary_data.append({
                    "File": filename, "Locus_ID": j+1,
                    "D0_Strength": round(final_popt[3*j], 4),
                    "r_Decay": round(final_popt[3*j+1], 3),
                    "Position": round(final_popt[3*j+2], 3),
                    "R_Squared": round(r2, 3)
                })
            print(f"成功处理: {filename} (R2={r2:.3f})")
        except Exception as e:
            print(f"处理失败 {file_path}: {e}")

    if summary_data:
        pd.DataFrame(summary_data).to_csv(os.path.join(args.outputdir, "all_chromosomes_detailed_loci.csv"), index=False)

    if all_results:
        num_files = len(all_results)
        fig, axes = plt.subplots(1, num_files, figsize=(6 * num_files, 6), sharey=True)
        if num_files == 1: axes = [axes]

        for i, (ax, (X, Y, popt, title, r2)) in enumerate(zip(axes, all_results)):
            ax.scatter(X, Y, color='gray', alpha=0.3, s=12, label="Observed")
            
            X_smooth = np.linspace(X.min(), X.max(), 1000)
            k = len(popt) // 3
            colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, k)) if k > 0 else []
            
            for j in range(k):
                d0, r, pos = popt[3*j : 3*j+3]
                y_comp = single_locus_model(X_smooth, d0, r, pos)
                
                # --- 核心修改部分：绘制位置引导线 ---
                ax.axvline(x=pos, color=colors[j], linestyle=':', lw=1.5, alpha=0.6)
                
                ax.plot(X_smooth, y_comp, '--', lw=1.2, color=colors[j], alpha=0.7)
                ax.fill_between(X_smooth, 0, y_comp, color=colors[j], alpha=0.05)
                
                # 调整标注位置，避免重合
                y_text = d0 + (0.05 if d0 > 0 else -0.15)
                ax.text(pos, y_text, f"L{j+1}\n{d0:.2f}", 
                        color=colors[j], fontsize=8, ha='center', fontweight='bold',
                        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=1))
            
            if k > 0:
                y_total = additive_model_k_loci(X_smooth, *popt)
                ax.plot(X_smooth, y_total, 'r-', lw=2.5, label="Global Fit")
            
            ax.axhline(0, color='black', lw=0.8)
            clean_title = title.split('.')[-2] if '.' in title else title
            ax.set_title(rf"{clean_title} ($R^2={r2:.3f}$)", fontweight='bold', fontsize=12)
            ax.set_xlabel("Position (cM/Mb)")
            if i == 0: ax.set_ylabel("Distortion Index (D)")
            ax.grid(True, alpha=0.15)
            ax.set_ylim(-1.1, 1.1)

        plt.tight_layout()
        plt.savefig(os.path.join(args.outputdir, "combined_sd_landscape.png"), dpi=200)
        plt.savefig(os.path.join(args.outputdir, "combined_sd_landscape.pdf"))
        print(f"\n[完成] 结果保存在文件夹: {args.outputdir}")

if __name__ == "__main__":
    main()
```

---

### Stage 2: Global Optimization via Levenberg–Marquardt Algorithm

Once seed parameters are obtained from Stage 1, all *k* loci are refined simultaneously using `scipy.optimize.curve_fit` with the multi-locus tanh-saturated model. This global optimization accounts for inter-locus interference during parameter estimation.

#### Boundary Constraints (Biological Plausibility)

| Parameter | Lower bound | Upper bound | Rationale |
|-----------|------------|-------------|-----------|
| *D₀* (selective strength) | −1.2 | 1.2 | Distortion index cannot exceed ±1; slight margin for numerical flexibility |
| *r* (decay rate) | 0.01 | 8.0 | Prevents near-zero (no-decay) and unrealistically fast decay |
| POS (position) | min(*X*) | max(*X*) | Constrain loci within the marker range |

