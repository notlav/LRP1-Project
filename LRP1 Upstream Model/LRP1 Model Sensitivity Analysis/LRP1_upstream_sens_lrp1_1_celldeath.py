# %% Sensitivity analysis weight of LRP1ag = 0.05
# Comparing high vs. low LRP1ag activity 

# %%
from math import isclose
import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
from IPython.display import display
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# Replace with Netflux files 
import LRP1_upstream as file
import LRP1_upstream_params as params

[speciesNames, tau, ymax, y0, w, n, EC50] = params.loadParams()

#%% Time span
tspan = [0, 50]

# Run default simulation to baseline
sol0 = solve_ivp(file.ODEfunc, tspan,  y0, args=(tau, ymax, w, n, EC50), t_eval=np.linspace(*tspan, 201))
sol0_SS = sol0.y[:,-1]

# Run simulation with perturbations to baseline
w_lrp1ag = 1
w_ros = 9

w[w_lrp1ag] = 1
# simulate MI 
w[w_ros] = 0.4

sol1 = solve_ivp(file.ODEfunc, tspan,  sol0_SS, args=(tau, ymax, w, n, EC50), t_eval=np.linspace(*tspan, 201))
sol1_SS = sol1.y[:,-1]

#%% ODE solutions into dataframe
sol0_df = pd.DataFrame(sol0.y.T, index = sol0.t, columns = speciesNames).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')
sol1_df = pd.DataFrame(sol1.y.T, index = sol1.t, columns = speciesNames).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')
                       
# create new column to label condition
sol0_df['condition'] = 'base'
sol1_df['condition'] = 'ROS, high lrp1ag'

# Here are the activity levels for each of your species for each condition
results = pd.concat([sol0_df, sol1_df])

# steady state values for each species
results_ss = results.query('time == @tspan[1]').reset_index().drop(columns=['time', 'index'])
display(results_ss)

#%% knockdown species at baseline and species at perturbation
sol_kd0 = np.zeros([len(ymax), len(ymax)]) 
sol_kd1 = np.zeros([len(ymax), len(ymax)])
sens = np.zeros([len(ymax), len(ymax)])

deltaP = -1 # complete knockdown

# simulating knockdowns of each species, one by one
for i in range(0, len(ymax)):
    print('Knockdown #', str(i+1), 'of', str(len(ymax)))
    ymax_kd = ymax.copy()
    # ymax_new = new maximum values of each species AFTER KNOCKDOWN
    # ymax_new should all be 0's since we are doing a complete knockdown
    ymax_kd[i] = (1+deltaP)*ymax[i]

    # baseline 
    sol0_kd = solve_ivp(file.ODEfunc, tspan, sol0_SS, args=(tau, ymax_kd, w, n, EC50))
    sol0_kd_SS = sol0_kd.y[:,-1]

    sol1_kd = solve_ivp(file.ODEfunc, tspan, sol0_kd_SS, args=(tau, ymax_kd, w, n, EC50))
    sol1_kd_SS = sol1_kd.y[:,-1]

    # change in activity for each species (KD at steady state - baseline at steady state)
    # Normalized sensitivity
    #sens[:,i] = -(sol1_kd_SS-sol0_SS)/(ymax_kd[i]-ymax[i])*ymax[i]/(sol0_SS)
    sens[:,i] = sol1_kd_SS-sol1_SS


# if there are any non-real numbers, change them to 0 
for i in range(0,len(ymax)):
    for j in range(0, len(ymax)):
        if np.isnan(np.real(sens[i,j]==1)):
            sens[i,j] == 0

#%% Plot Sensitivity Matrix for all species

vcenter = 0
vmin, vmax = -1, 1 #sens_kd.min(), sens_kd.max()
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
#colormap = cm.RdBu_r
colormap = 'seismic'

sens_heatmap = pd.DataFrame(sens)
fig1, ax1 = plt.subplots()
fig1.set_size_inches(20.5, 14.5)
ax1 = sns.heatmap(sens_heatmap, norm=normalize, cmap=colormap, xticklabels=speciesNames, yticklabels=speciesNames)

font = {'fontname':'Arial'}
#plt.title("Sensitivity Analysis", fontsize = 40, **font)
plt.xlabel("Node Knockdown", fontsize = 30, **font)
plt.ylabel("Measured Node Activity", fontsize=30, **font)

#%% export sensitivity matrix
sens_df = sens_heatmap.set_index([speciesNames]).set_axis([speciesNames], axis = 1)
#sens_df.to_csv('/Volumes/SaucermanLab/Lavie/LRP1/Code/LRP1 Upstream Model/sens_lrp1_0.csv')

# %% Find indices of species
speciesNames.index('cellDeath')

#%% plot cell death sensitivity only 

sens_celldeath = pd.DataFrame(sens[12,:])
sens_celldeath.index = speciesNames
sens_celldeath.columns = ['cellDeath activity']

fig1, ax1 = plt.subplots()
fig1.set_size_inches(8,6)

# Filter species that have cellDeath activity > 0.1
celldeath_nodes = sens_celldeath[sens_celldeath['cellDeath activity'] > 0.1]

# sort ranked nodes
ranked_nodes = celldeath_nodes.sort_values(by='cellDeath activity', ascending=False)

ax1 = sns.barplot(ranked_nodes.transpose(), color = 'black')
# rotate x labels
plt.xticks(rotation=45)

font = {'fontname':'Arial'}
#plt.title("Sensitivity Analysis", fontsize = 40, **font)
plt.xlabel("Node Knockdown", fontsize = 20, **font)
plt.ylabel("Cell Death (Change in Activity)", fontsize=20, **font)

# export as svg
plt.savefig('/Volumes/SaucermanLab/Lavie/LRP1/Figures/sensCellDeath.svg', format='svg', dpi=1200)

plt.show()

# %%
