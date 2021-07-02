import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 

# Figure 6B
# Relative frequency of each variant of each sub-library
# Python script was run in JupyterLab

# set file name
save = 'fig_unique_freq.png'

# set p to 100 to show y-axis as percent
# set p to 1 to show y-axis as a decimal
p = 1

# Read pkl file
pklfile = 'mdf_FN80cysaCA2.pkl'
ot = pd.read_pickle(pklfile)
pklfile = 'mdf_FN28cysaCA2.pkl'
ot2 = pd.read_pickle(pklfile)
pklfile = 'mdf_FN28cysaCA9.pkl'
ot3 = pd.read_pickle(pklfile)

# Drop counts for naive libraries
ot = ot.drop(['C80_0'], axis=1)  #Fn80
ot2 = ot2.drop(['C28_0'], axis=1) #Fn28
ot3 = ot3.drop(['C28_0'], axis=1) #Fn28aCA9

# Only include sequences that are present in the sub-library
Fn80p2 = ot[ot['C80P2']>0].nlargest(494, 'C80P2')['C80P2']
Fn80p3 = ot[ot['C80P3']>0].nlargest(205, 'C80P3')['C80P3']
Fn80p5 = ot[ot['C80P5']>0].nlargest(169, 'C80P5')['C80P5']
Fn80p7 = ot[ot['C80P7']>0].nlargest(180, 'C80P7')['C80P7']

Fn28p2 = ot2[ot2['C28P2CA2']>0].nlargest(199, 'C28P2CA2')['C28P2CA2']
Fn28p3 = ot2[ot2['C28P3CA2']>0].nlargest(12, 'C28P3CA2')['C28P3CA2']
Fn28p5 = ot2[ot2['C28P5CA2']>0].nlargest(429, 'C28P5CA2')['C28P5CA2']
Fn28p7 = ot2[ot2['C28P7CA2']>0].nlargest(443, 'C28P7CA2')['C28P7CA2']

Fn28p2aCA9 = ot3[ot3['C28P2CA9']>0].nlargest(47, 'C28P2CA9')['C28P2CA9']
Fn28p3aCA9 = ot3[ot3['C28P3CA9']>0].nlargest(83, 'C28P3CA9')['C28P3CA9']
Fn28p5aCA9 = ot3[ot3['C28P5CA9']>0].nlargest(90, 'C28P5CA9')['C28P5CA9']
Fn28p7aCA9 = ot3[ot3['C28P7CA9']>0].nlargest(421, 'C28P7CA9')['C28P7CA9']

# Determine relative frequency value
Fn80p2 = (Fn80p2.div(Fn80p2.sum())).to_numpy() * p
Fn80p3 = (Fn80p3.div(Fn80p3.sum())).to_numpy() * p
Fn80p5 = (Fn80p5.div(Fn80p5.sum())).to_numpy() * p
Fn80p7 = (Fn80p7.div(Fn80p7.sum())).to_numpy() * p
Fn28p2 = (Fn28p2.div(Fn28p2.sum())).to_numpy() * p
Fn28p3 = (Fn28p3.div(Fn28p3.sum())).to_numpy() * p
Fn28p5 = (Fn28p5.div(Fn28p5.sum())).to_numpy() * p
Fn28p7 = (Fn28p7.div(Fn28p7.sum())).to_numpy() * p
Fn28p2aCA9 = (Fn28p2aCA9.div(Fn28p2aCA9.sum())).to_numpy() * p
Fn28p3aCA9 = (Fn28p3aCA9.div(Fn28p3aCA9.sum())).to_numpy() * p
Fn28p5aCA9 = (Fn28p5aCA9.div(Fn28p5aCA9.sum())).to_numpy() * p
Fn28p7aCA9 = (Fn28p7aCA9.div(Fn28p7aCA9.sum())).to_numpy() * p

# Values for x-axis, # of unique variants
Fn80p2x = range(1, len(Fn80p2)+1, 1)
Fn80p3x = range(1, len(Fn80p3)+1, 1)
Fn80p5x = range(1, len(Fn80p5)+1, 1)
Fn80p7x = range(1, len(Fn80p7)+1, 1)
Fn28p2x = range(1, len(Fn28p2)+1, 1)
Fn28p3x = range(1, len(Fn28p3)+1, 1)
Fn28p5x = range(1, len(Fn28p5)+1, 1)
Fn28p7x = range(1, len(Fn28p7)+1, 1)
Fn28p2aCA9x = range(1, len(Fn28p2aCA9)+1, 1)
Fn28p3aCA9x = range(1, len(Fn28p3aCA9)+1, 1)
Fn28p5aCA9x = range(1, len(Fn28p5aCA9)+1, 1)
Fn28p7aCA9x = range(1, len(Fn28p7aCA9)+1, 1)


# Generation of figure
fig, ax = plt.subplots(figsize=(5,5), dpi=150)

plt.plot(Fn28p2aCA9x,Fn28p2aCA9, color='#1f77b4', linewidth = 4, linestyle = 'solid')
plt.plot(Fn28p3aCA9x,Fn28p3aCA9, color='#ff7f0e', linewidth = 4, linestyle = 'solid') 
plt.plot(Fn28p5aCA9x,Fn28p5aCA9, color='#2ca02c', linewidth = 4, linestyle = 'solid') 
plt.plot(Fn28p7aCA9x,Fn28p7aCA9, color='#d62728', linewidth = 4, linestyle = 'solid') 

plt.plot(Fn28p2x,Fn28p2, color='#1f77b4', linewidth = 2, linestyle = 'solid') 
plt.plot(Fn28p3x,Fn28p3, color='#ff7f0e', linewidth = 2, linestyle = 'solid') 
plt.plot(Fn28p5x,Fn28p5, color='#2ca02c', linewidth = 2, linestyle = 'solid') 
plt.plot(Fn28p7x,Fn28p7, color='#d62728', linewidth = 2, linestyle = 'solid')

plt.plot(Fn80p2x,Fn80p2, color='#1f77b4', linewidth = 2, linestyle = 'dashed') 
plt.plot(Fn80p3x,Fn80p3, color='#ff7f0e', linewidth = 2, linestyle = 'dashed') 
plt.plot(Fn80p5x,Fn80p5, color='#2ca02c', linewidth = 2, linestyle = 'dashed') 
plt.plot(Fn80p7x,Fn80p7, color='#d62728', linewidth = 2, linestyle = 'dashed') 

#ax.legend(['PEG2','PEG3','PEG5','PEG7'], fontsize=10, loc = 'right')

plt.xlim([1, 1000])
plt.ylim([0.000001, 1])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Unique Variant', fontsize=16)
plt.ylabel('% of Reads', fontsize=16)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

plt.show()
fig.savefig(save, bbox_inches = "tight")
