import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sim_var_same = pd.read_csv('sim_df/variation_umi_stripplot.tsv', sep='\t')
sim_var_80eff = pd.read_csv('sim_df/umi_stripplot_80eff.tsv', sep='\t')
sim_var_90eff = pd.read_csv('sim_df/umi_stripplot_90eff.tsv', sep='\t')
sim_var_95eff = pd.read_csv('sim_df/umi_stripplot_95eff.tsv', sep='\t')

sim_var_5cycles = pd.read_csv('sim_df/umi_stripplot_5cycles.tsv', sep='\t')
sim_var_10cycles = pd.read_csv('sim_df/umi_stripplot_10cycles.tsv', sep='\t')
sim_var_20cycles = pd.read_csv('sim_df/umi_stripplot_20cycles.tsv', sep='\t')

plt.figure(figsize=(5, 5))
sns.swarmplot(data=sim_var_same, x='Simulation', y='Occurrences',color='.2')
plt.xticks([])
plt.title('Violin Plot of String Occurrences (Subset of 30 Categories)')
plt.savefig('variation_umi_stripplot.png')

