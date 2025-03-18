import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

# Reading data from the saved file
df = pd.read_csv("last_nt_time.csv")

# Preparing DataFrames with a copy to avoid warnings
df_last_time = df[['Mzams', 'E_final', 'last_time']]
df_last_time['last_time'] = df_last_time['last_time'] / 86400  # Convert seconds to days
df_last_time_pivot = df_last_time.pivot('Mzams', 'E_final', 'last_time')

df_last_nt = df[['Mzams', 'E_final', 'last_nt']]
df_last_nt_pivot = df_last_nt.pivot('Mzams', 'E_final', 'last_nt')


# Function for custom palette
def mix_palette():
    palette = sns.color_palette("viridis", 30)
    palette[0] = 'indigo'
    return palette


# Plotting heatmaps (steps)
plt.figure(figsize=(8, 6))
ax1 = sns.heatmap(
    df_last_nt_pivot,
    cmap=mix_palette(),
    square=True,
    xticklabels=True,
    yticklabels=True,
    cbar_kws={'label': 'Steps'}  # Label for colorbar
)
ax1.figure.axes[-1].yaxis.label.set_size(26)  # Adjust font size for colorbar label
ax1.figure.axes[-1].tick_params(labelsize=18)  # Adjust font size for colorbar ticks
colorbar = ax1.collections[0].colorbar
colorbar.ax.yaxis.offsetText.set_fontsize(16)  # Adjust font size for the scientific notation label
plt.xlabel(r"$E_{\mathrm{final}}$", fontsize=36)
plt.ylabel(r"$M_{\mathrm{ZAMS}}$", fontsize=36)
plt.xticks(fontsize=18)  # Adjust font size for x-axis ticks
plt.yticks(fontsize=18)  # Adjust font size for y-axis ticks
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()


# Plotting heatmaps (days)
plt.figure(figsize=(8, 6))
ax2 = sns.heatmap(
    df_last_time_pivot,
    cbar_kws={'label': 'Days'}  # Label for colorbar
)
ax2.figure.axes[-1].yaxis.label.set_size(26)  # Adjust font size for colorbar label
ax2.figure.axes[-1].tick_params(labelsize=18)  # Adjust font size for colorbar ticks
plt.xlabel(r"$E_{\mathrm{final}}$", fontsize=36)
plt.ylabel(r"$M_{\mathrm{ZAMS}}$", fontsize=36)
plt.xticks(fontsize=18)  # Adjust font size for x-axis ticks
plt.yticks(fontsize=18)  # Adjust font size for y-axis ticks
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.show()