# Create bar plots from validation results
# author: Lavie
# date: 7/1/2024

#%% import statements
import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import LRP1_upstream as LRP1_ODE
import LRP1_upstream_params as LRP1_ODE_params

#%% run model simulation

# LRP1 low and high conditions
LRP1_low = 0
LRP1_high = 1

# reaction indices for LRP1 and ROS
idx_LRP1ag = 0
idx_ROS = 9

[speciesNames, tau, ymax, y0, w, n, EC50] = LRP1_ODE_params.loadParams()

tspan = [0, 50]

# simulate MI
w[idx_ROS] = 0.4

# simulate low LRP1ag
w[idx_LRP1ag] = LRP1_low
sol_low= solve_ivp(LRP1_ODE.ODEfunc, tspan, y0, args=(tau, ymax, w, n, EC50), t_eval=np.linspace(*tspan, 201))
lowLRP1_SS = sol_low.y[:,-1]

# simulate high LRP1
w[idx_LRP1ag] = LRP1_high
sol_high = solve_ivp(LRP1_ODE.ODEfunc, tspan, lowLRP1_SS, args=(tau, ymax, w, n, EC50), t_eval=np.linspace(*tspan, 201))
highLRP1_SS = sol_high.y[:,-1]

#%% create dataframes for steady state values for each species
lowLRP1_df = pd.DataFrame(sol_low.y.T, index = sol_low.t, columns = speciesNames).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')

highLRP1_df = pd.DataFrame(sol_high.y.T, index = sol_high.t, columns = speciesNames).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')

# grab steady state values for each species (time = 50)
# filter time = 50
lowLRP1_ss = lowLRP1_df.query('time == @tspan[1]').reset_index().drop(columns=['time', 'index'])

highLRP1_ss = highLRP1_df.query('time == @tspan[1]').reset_index().drop(columns=['time', 'index'])

# combine low and high LRP1 steady state values and add condition column
lowLRP1_ss['condition'] = 'control'
highLRP1_ss['condition'] = 'LRP1\nagonist'

steady_long = pd.concat([lowLRP1_ss, highLRP1_ss])


#%% plot activity for each output
# plot all species in one figure 
validation_species = ('ERK12', 'PI3K', 'Akt', 'NFkB', 'Bax', 'Bcl2', 'cas3', 'apoptosis')

fig, ax = plt.subplots(2,4, figsize = (20, 10), sharey = True)

palette = sns.color_palette('crest', n_colors=8)

# plot each species in each subplot
for i, species in enumerate(validation_species):
	sns.barplot(
		data = steady_long[steady_long['species'] == species],
		x = 'condition',
		y = 'activity',
		color= palette[i],
		ax = ax[i//4, i%4]
	)
	ax[i//4, i%4].set_title(f'{species} activity', fontsize = 16)
	ax[i//4, i%4].set_xlabel('')
	ax[i//4, i%4].set_ylabel('activity')
	# set font to arial
	plt.rcParams['font.family'] = 'Arial'

#plt.savefig('/Volumes/SaucermanLab/Lavie/LRP1/Figures/' + 'validation_activity_gradient.svg')

#%% plot activity for each output under different ROS conditions
fig, axs = plt.subplots()

# loop through each species and plot activity for each condition
# show each plot separately
# add title
palette = sns.color_palette('crest', n_colors=8)

for i, species in enumerate(validation_species):
	plt.figure(figsize = (3, 3))

	sns.barplot(
		data = steady_long[steady_long['species'] == species],
		x = 'condition',
		y = 'activity',
		color= palette[i],
		width = 0.5,
	)

	plt.xticks(fontsize = 14)
	plt.yticks(fontsize = 14)
	plt.ylabel('activity', fontsize = 16)
	plt.xlabel('', fontsize = 16)

	plt.title(f'{species} activity', fontsize = 16)
	# set font to arial
	plt.rcParams['font.family'] = 'Arial'

	# export to svg
	# change file name each time
	# include file path
	#plt.savefig('/Volumes/SaucermanLab/Lavie/LRP1/Figures/' + f'{species}_activity.svg')

	plt.show()
# %% manually plot experimental cas3 data
exp_cas3 = pd.DataFrame({'condition': ['control', 'LRP1\nagonist'], 'activity': [2.05, 0.57]})
model_cas3 = steady_long[steady_long['species'] == 'cas3']

# plot subplot for exp and model cas3
fig, axes = plt.subplots(1, 2, figsize = (8, 3),)
plt.suptitle('caspase 3 activity', fontsize = 16)

color_palette = sns.color_palette('crest', n_colors=8)
# model prediction
sns.barplot(model_cas3, x = 'condition', y = 'activity', width = 0.5, color = color_palette[6], ax = axes[0])
axes[0].set_xlabel('Model', fontsize = 16)
axes[0].set_ylabel('activity', fontsize = 12)
axes[0].set_yticks(np.arange(0, 1.2, 0.2))
axes[0].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

# experimental cas3 data
sns.barplot(exp_cas3, x = 'condition', y = 'activity', yerr = [0.495, 0.4], width = 0.5, color = color_palette[6], ax = axes[1])
axes[1].set_xlabel('Experiment', fontsize = 16)
axes[1].set_ylabel('-fold increase (vs. sham)', fontsize = 12)
axes[1].set_yticks(np.arange(0, 3, 0.5))
axes[1].set_yticklabels([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

# set font name to arial
plt.show()

plt.savefig('Y:/Lavie/LRP1/Figures/' + 'cas3ValidationBarplot_final.svg', bbox_inches = 'tight')


# %%
