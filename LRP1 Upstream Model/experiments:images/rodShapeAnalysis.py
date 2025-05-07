# author: LLN
# date: 8/10/24
# create barplot for % rod shape analysis generated from cellProfiler

#%% import libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

#%% import excel data
shapeData = pd.read_excel('/Volumes/SaucermanLab/Lavie/LRP1/Experiments/CellProfiler/PercentRodOutputs/RodShapeAnalysis.xlsx')

#%% calculate mean and std for each condition
conditions = ['Normoxia', 'I-R', 'I-R \n+ \nDrug 1', 'I-R \n+ \nDrug 2']

# assign data to each condition
normoxia = shapeData[shapeData['treatment'].str.contains('normoxia', case = False)]
hypoxia = shapeData[shapeData['treatment'].str.contains('hypoxia', case = False)]
SP16 = shapeData[shapeData['treatment'].str.contains('SP16', case = False)]
A2M = shapeData[shapeData['treatment'].str.contains('A2M', case = False)]

# calculate mean and std for each condition in new dataframe
avgs = [normoxia['% rod'].mean()*100, hypoxia['% rod'].mean()*100, SP16['% rod'].mean()*100, A2M['% rod'].mean()*100] 

stds = [normoxia['% rod'].std()*100, hypoxia['% rod'].std()*100, SP16['% rod'].std()*100, A2M['% rod'].std()*100]

#%% generate bar plot
plt.figure(figsize = (4, 2))

color_palette = sns.color_palette('crest', n_colors=8)

plt.bar(x = conditions, height = avgs, yerr = stds, color = color_palette[6], width = 0.4)

plt.ylabel('% Rod Shape', fontsize = 12)
plt.rcParams['font.family'] = 'Times New Roman'

# increase border height
plt.ylim(0, 100)
plt.show()

#%% run anova
f_val, p_val = stats.f_oneway(hypoxia['% rod'], normoxia['% rod'], SP16['% rod'], A2M['% rod'])

# tukey's post hoc test
data = np.concatenate([hypoxia['% rod'], normoxia['% rod'],  SP16['% rod'], A2M['% rod']])
labels = np.concatenate([['hypoxia']*len(hypoxia), ['normoxia']*len(normoxia), ['SP16']*len(SP16), ['A2M']*len(A2M)])

pt = pairwise_tukeyhsd(data, labels, alpha = 0.05)
print(pt)


# %%
