# Converted from 3_L2FC_visualisation.ipynb

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set the working directory
working_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords"
os.chdir(working_dir)
print(f"Current working directory: {os.getcwd()}")

# Load the Peak.csv file
peak_df = pd.read_csv('./Peak.csv')

# Display the first few rows to verify the data
# Remove the first row from the dataframe
peak_df = peak_df.iloc[1:]

# Display the first few rows to verify the data
peak_df.head()

# Remove all rows with negative L2FC
peak_df = peak_df[peak_df['L2FC'] >= 0]

print(peak_df.shape)

selected = ["GEMIN4", "YME1L1", "VCP", "DARS1", "POLR2A", "ORAI1", "POLRMT", "CBX2", "NSD2", "TLE1", "CACNG1", "ZMYND11", "DCAF5", "MED25", "SEMA6A-AS2", "POLD3", "RFT1", "SETDB1", "POLE", "CHD4", "POLR2C", "HTT", "MED14", "PER3", "TOR3A", "EEF1A2", "IL33", "SETD5", "GNRH1", "POLR2E", "IGF1R", "IGF2R", "MTMR1", "EEF2", "RBPMS", "EIF2S3B", "NIPBL", "VCP"]



# Sort the DataFrame by L2FC values
sorted_df = peak_df.sort_values('L2FC', ascending=False)

# Select top 20 genes by absolute L2FC for better readability in some plots
sorted_df_20 = sorted_df.iloc[:20]

# Select top 50 genes by absolute L2FC for better readability in some plots
sorted_df_50 = sorted_df.iloc[:50]

# Select top 150 genes by absolute L2FC for better readability in some plots
sorted_df_150 = sorted_df.iloc[:150]


plt.figure(figsize=(12, 6))
bars = sns.barplot(x='SYMBOL', y='L2FC', data=sorted_df_50)
plt.title('Top 50 Genes by Log2 Fold Change')
plt.xticks(rotation=45, ha='right')

# Color bars in red if their symbol is in the selected list
for i, bar in enumerate(bars.patches):
    if sorted_df_50.iloc[i]['SYMBOL'] in selected:
        bar.set_facecolor('red')

plt.tight_layout()
plt.show()

# Create a figure with subplots (6 rows, 4 columns)
fig, axes = plt.subplots(6, 4, figsize=(24, 36))
fig.suptitle('Top 50 Genes by Log2 Fold Change for Each Chromosome', fontsize=16, y=1.02)

# Flatten the axes array for easier iteration
axes = axes.flatten()

chromosomes = ['chr2', 'chr17', 'chr12', 'chr10', 'chr7', 'chr9', 'chr19',
       'chr13', 'chr22', 'chr6', 'chr1', 'chr4', 'chr11', 'chr5', 'chr8',
       'chr3', 'chr14', 'chr16', 'chr18', 'chrX', 'chr15', 'chr21', 'chr20']

for idx, seqname in enumerate(chromosomes):
    # Filter data for the current seqname
    seqname_df = sorted_df[sorted_df['seqnames'] == seqname].head(50)
    
    # Create subplot
    bars = sns.barplot(x='SYMBOL', y='L2FC', data=seqname_df, ax=axes[idx])
    
    # Color bars in red if their symbol is in the selected list
    for i, bar in enumerate(bars.patches):
        if seqname_df.iloc[i]['SYMBOL'] in selected:
            bar.set_facecolor('red')
    
    axes[idx].set_title(f'{seqname}')
    axes[idx].set_xticklabels([])  # Remove x-axis labels for better readability
    axes[idx].set_xlabel('')  # Remove x-axis label
    
    if idx % 4 == 0:  # Only set y-label for leftmost subplots
        axes[idx].set_ylabel('Log2 Fold Change')
    else:
        axes[idx].set_ylabel('')

# Remove any unused subplots
for idx in range(len(chromosomes), 24):
    fig.delaxes(axes[idx])

plt.tight_layout()
plt.show()

def plot_l2fc_distribution(sorted_df, n_top=5, selected=None):
    plt.figure(figsize=(16, 8)) 

    # Plot all points
    sns.scatterplot(x=range(len(sorted_df)), y='L2FC', data=sorted_df, color='gray', alpha=0.5)

    # Highlight selected points in green
    if selected is not None:
        selected_df = sorted_df[sorted_df['SYMBOL'].isin(selected)]
        sns.scatterplot(x=selected_df.index, y='L2FC', data=selected_df, color='green', s=50)

    # Highlight top n_top points
    top_n = sorted_df.head(n_top)
    sns.scatterplot(x=top_n.index, y='L2FC', data=top_n, color='red', s=100)

    # Add labels for top n_top points with improved positioning
    label_offset_y = 0.2  # Vertical offset between labels
    label_offset_x = 50   # Horizontal offset for labels

    for i, (idx, row) in enumerate(top_n.iterrows()):
        x = idx
        y = row['L2FC']
        label = row['SYMBOL']
        
        # Calculate vertical and horizontal position for the label
        y_offset = label_offset_y * (n_top - i)  # (i - n_top // 2) * label_offset_y
        if y_offset > 0.4:
            y_offset = 0.4
        x_offset = label_offset_x * (-1 if i % 2 == 0 else 1)  # Alternate left and right
        
        plt.annotate(label, (x, y), 
                     xytext=(x + x_offset, y + y_offset), 
                     textcoords='data',
                     ha='left' if x_offset > 0 else 'right', va='center',
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1'))

    plt.title(f'Log2 Fold Change Distribution with Top {n_top} Highlighted', fontsize=16)
    plt.xlabel('Gene Index', fontsize=12)
    plt.ylabel('Log2 Fold Change', fontsize=12)

    # Ensure L2FC values are numeric and remove any non-numeric values
    sorted_df['L2FC'] = pd.to_numeric(sorted_df['L2FC'], errors='coerce')
    min_l2fc = np.floor(sorted_df['L2FC'].min())
    max_l2fc = np.ceil(sorted_df['L2FC'].max())

    # Calculate appropriate y-ticks
    y_ticks = np.linspace(min_l2fc, max_l2fc, 6)
    plt.yticks(y_ticks)

    # Improve x-axis
    plt.xlim(-50, len(sorted_df) + 50)  # Add some padding
    x_ticks = np.linspace(0, len(sorted_df), 5, dtype=int)
    plt.xticks(x_ticks, x_ticks)

    plt.tick_params(axis='both', which='major', labelsize=10)

    plt.tight_layout()
    plt.show()

plot_l2fc_distribution(sorted_df, n_top=5, selected=selected)

from matplotlib.lines import Line2D

plt.figure(figsize=(12, 6))

# Create the main histogram
sns.histplot(data=sorted_df, x='L2FC', kde=True, color='lightblue', edgecolor='black')

# Highlight the selected symbols
for symbol in selected:
    l2fc = sorted_df[sorted_df['SYMBOL'] == symbol]['L2FC'].values[0]
    plt.axvline(x=l2fc, color='red', linestyle='--', alpha=0.7)
    # plt.text(l2fc, plt.gca().get_ylim()[1], symbol, rotation=90, va='top', ha='right', color='red')

plt.title('Distribution of Log2 Fold Change with Selected Symbols Highlighted')
plt.xlabel('Log2 Fold Change')
plt.ylabel('Count')

# Add a legend
custom_lines = [Line2D([0], [0], color='lightblue', lw=4),
                Line2D([0], [0], color='red', linestyle='--', lw=2)]
plt.legend(custom_lines, ['All Genes', 'Selected Symbols'])

plt.tight_layout()
plt.show()



