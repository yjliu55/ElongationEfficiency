import os
import glob
import re
import numpy as np
import pandas as pd
import pyBigWig
from pyfaidx import Fasta
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from adjustText import adjust_text
except ImportError:
    print("Please install adjustText: pip install adjustText")
    adjust_text = None

# ================= 配置路径 =================
TXT_FILE = "merge_WT_22_rep_average.txt"
FASTA_FILE = "Arabidopsis.fa"
BW_DIR = "./"  

# ================= 1. 读取目标数据 =================
print("Loading target data...")
df = pd.read_csv(TXT_FILE, sep=r'\s+')
df.rename(columns={'0': 'elongation_efficiency'}, inplace=True)
df['chrom'] = df['chrom'].astype(str) 

median_eff = df['elongation_efficiency'].median()
df['efficiency_class'] = (df['elongation_efficiency'] >= median_eff).astype(int)

# ================= 2. 提取序列特征 (新增 AT 含量) =================
print("Extracting sequence contexts...")
fasta = Fasta(FASTA_FILE)

def get_motif_frequencies(chrom, start, end):
    try:
        seq = fasta[chrom][start:end].seq.upper()
        L = len(seq)
        if L == 0:
            return pd.Series({m: 0.0 for m in ['CG', 'CHG', 'CHH', 'GC', 'AT']})
        
        counts = {
            'CG': len(re.findall(r'(?=CG)', seq)),
            'CHG': len(re.findall(r'(?=C[ACT]G)', seq)),
            'CHH': len(re.findall(r'(?=C[ACT][ACT])', seq)),
            'GC': seq.count('G') + seq.count('C'),
            'AT': seq.count('A') + seq.count('T')  # 新增 AT 含量特征
        }
        freqs = {k: v / L for k, v in counts.items()}
        return pd.Series(freqs)
    except KeyError:
        return pd.Series({m: np.nan for m in ['CG', 'CHG', 'CHH', 'GC', 'AT']})

motif_df = df.apply(lambda row: get_motif_frequencies(row['chrom'], row['start'], row['end']), axis=1)
df = pd.concat([df, motif_df], axis=1)

# ================= 3. 提取表观修饰信号 =================
print("Extracting epigenetic signals...")
bw_files = glob.glob(os.path.join(BW_DIR, "*.bw"))
bw_features = []

for bw_path in bw_files:
    bw_name = os.path.basename(bw_path).replace('.bw', '').replace('_sort', '')
    bw_features.append(bw_name)
    try:
        bw = pyBigWig.open(bw_path)
        def get_mean_signal(chrom, start, end):
            try:
                val = bw.stats(chrom, start, end, type="mean")[0]
                return val if val is not None else 0.0
            except RuntimeError:
                return 0.0
        df[bw_name] = df.apply(lambda row: get_mean_signal(row['chrom'], row['start'], row['end']), axis=1)
        bw.close()
    except Exception as e:
        print(f"Error reading {bw_path}: {e}")

# ================= 4. 特征分类 =================
feature_categories = {}
sequence_features = ['CG', 'CHG', 'CHH', 'GC', 'AT', 'effective_length']
for f in sequence_features:
    feature_categories[f] = 'others'
for f in bw_features:
    if 'CG' in f or 'CHH' in f or 'CHG' in f or 'methylation' in f.lower():
        feature_categories[f] = 'DNA methylation'
    else:
        feature_categories[f] = 'histone modifications & H1'

# ================= 5. 数据清洗与建模 =================
print("Preparing model...")
df = df.dropna()
X = df[sequence_features + bw_features]
y = df['efficiency_class']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

rf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)

# ================= 6. 计算 Importance =================
print("Calculating Feature Importances...")
gini_decrease = rf.feature_importances_ * X_train.shape[0]

perm_imp = permutation_importance(rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=-1)
accuracy_decrease = perm_imp.importances_mean

# ================= 7. P-value 置换检验 (优化计算逻辑) =================
print("Running permutation test for p-values...")
n_permutations = 100
null_acc = np.zeros((n_permutations, X.shape[1]))

y_train_array = y_train.values.copy()
for i in range(n_permutations):
    np.random.shuffle(y_train_array)
    rf_null = RandomForestClassifier(n_estimators=50, random_state=i, n_jobs=-1)
    rf_null.fit(X_train, y_train_array)
    
    perm_null = permutation_importance(rf_null, X_test, y_test, n_repeats=1, random_state=i, n_jobs=-1)
    null_acc[i, :] = perm_null.importances_mean

# 引入标准的经验 p 值计算公式 (Count + 1) / (N + 1)
p_values = (np.sum(null_acc >= accuracy_decrease, axis=0) + 1) / (n_permutations + 1)

plot_df = pd.DataFrame({
    'Feature': X.columns,
    'accuracy_decrease': accuracy_decrease,
    'gini_decrease': gini_decrease,
    'p_value': p_values,
    'Category': [feature_categories[f] for f in X.columns]
})
plot_df['significance'] = np.where(plot_df['p_value'] < 0.01, 'p < 0.01', 'p >= 0.01')

# ================= 8. 绘制图表 (强制颜色映射) =================
print("Plotting results...")
plt.figure(figsize=(7, 7), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')
ax.grid(True, color='#e5e5e5', linestyle='-', linewidth=1)
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(1)

# 确保配色方案的绝对映射
colors_pval = {'p < 0.01': '#F8766D', 'p >= 0.01': '#00BFC4'}

# 如果本次运算有 < 0.01 的数据点，它们会被准确映射为红色
sns.scatterplot(
    data=plot_df, 
    x='accuracy_decrease', 
    y='gini_decrease', 
    hue='significance',
    hue_order=['p < 0.01', 'p >= 0.01'], # 强制图例和着色顺序
    palette=colors_pval,
    s=80, 
    edgecolor='white',
    linewidth=0.5,
    alpha=0.9,
    ax=ax
)

colors_text = {
    'DNA methylation': '#00BFFF',
    'histone modifications & H1': '#FF69B4',
    'others': '#000000'
}

texts = []
for idx, row in plot_df.iterrows():
    display_name = row['Feature']
    
    display_name = re.sub(r'_\d+$', '', display_name) 
    if display_name.startswith('WT_'):
        display_name = display_name[3:]
    if 'average' in row['Feature'] and ('CG' in row['Feature'] or 'CH' in row['Feature']):
         display_name = 'm' + display_name.replace('_average', '')
         
    t = ax.text(
        row['accuracy_decrease'], 
        row['gini_decrease'], 
        display_name,
        color=colors_text[row['Category']],
        fontsize=11,
        fontweight='medium'
    )
    texts.append(t)

if adjust_text:
    adjust_text(texts, 
                x=plot_df['accuracy_decrease'].values,
                y=plot_df['gini_decrease'].values,
                force_text=(1.0, 1.5), 
                force_points=(1.0, 1.5), 
                expand_points=(2.0, 2.0),
                max_iterations=1500,
                arrowprops=dict(arrowstyle="-", color='gray', lw=0.5, alpha=0.7))

ax.set_xlabel('accuracy decrease', fontsize=14, color='black')
ax.set_ylabel('gini decrease', fontsize=14, color='black')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, loc='upper left', frameon=False, fontsize=11, handletextpad=0.1)

plt.text(0.98, 1.05, 'DNA methylation', color=colors_text['DNA methylation'], 
         ha='right', va='top', transform=ax.transAxes, fontsize=11, fontweight='bold')
plt.text(0.98, 1.00, 'histone modifications & H1', color=colors_text['histone modifications & H1'], 
         ha='right', va='top', transform=ax.transAxes, fontsize=11, fontweight='bold')
plt.text(0.98, 0.95, 'others', color=colors_text['others'], 
         ha='right', va='top', transform=ax.transAxes, fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig('RF_Feature_Importance_Plot_Optimized.pdf', format='pdf', bbox_inches='tight')
print("Analysis complete. Plot saved as 'RF_Feature_Importance_Plot_Optimized.pdf'.")
