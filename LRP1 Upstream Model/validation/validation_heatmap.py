# Lavie Ngo
# Creates 2-column heatmap to compare input/output computational vs. experimental results from literature for validation
# Categories: Increase, Decrease, No change

#%% Import statements
import pandas as pd
import numpy as np
import openpyxl
import matplotlib.pyplot as plt 
import seaborn as sns
from IPython.display import display
import matplotlib.colors as mcolors
import matplotlib.cm as cm

#%% Read validation excel file
validation = pd.read_excel('/Volumes/SaucermanLab/Lavie/LRP1/Code/LRP1 Upstream Model/validation/LRP1_validation_edit.xlsx')

# %% Take 'Measurement' column and create a new column mapping to numerical values 
lit = np.zeros(len(validation))
model = np.zeros(len(validation))

# Literature results
for i in range(0,len(lit)):
    if validation['Measurement'][i] == 'Increase':
        lit[i] = 1
    elif validation['Measurement'][i] == 'Decrease':
        lit[i] = -1 
    elif validation['Measurement'][i] == 'No Change':
        lit[i] = 0

# Model results
for i in range(0,len(model)):
    if validation['Prediction'][i] == 'Increase':
        model[i] = 1
    elif validation['Prediction'][i] == 'Decrease':
        model[i] = -1 
    elif validation['Prediction'][i] == 'No Change':
        model[i] = 0

# add new column
validation['Literature Number Value'] = lit
validation['Model Number Value'] = model

results = {'Output': validation['Output'], 'Model': model, 'Experiments': lit}

validation_heatmap = pd.DataFrame(results).set_index('Output')
# sort 
validation_heatmap = validation_heatmap.sort_values('Model')

# select outputs
validation_heatmap_sub = validation_heatmap.loc[['Akt', 'cas3', 'apoptosis']]     
# %% Plot heatmap

vcenter = 0
vmin, vmax = -1, 1 
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

#colors = ["#497CB7", "#A42115"]
#colormap = sns.color_palette(colors)

fig1, ax1 = plt.subplots()
fig1.set_size_inches(8, 20)
ax1 = sns.heatmap((validation_heatmap), linewidths=1, norm=normalize, cmap='seismic')


font = {'fontname':'Arial'}
#plt.title("Sensitivity Analysis", fontsize = 40, **font)
plt.xlabel("", fontsize = 30, **font)
plt.ylabel("", fontsize=30, **font)
plt.yticks(rotation = 0, fontsize = 35, **font)
plt.xticks(fontsize = 35, **font)


# %%
