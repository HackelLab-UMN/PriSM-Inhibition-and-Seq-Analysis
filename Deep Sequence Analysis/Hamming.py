import numpy as np
import pandas as pd
import scipy.spatial.distance as sci
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.ticker import FormatStrFormatter

# Figure 6A
# Determination of Hamming distance
# Python script was used in JupyterLab

# Save figures as . . . 
save = 'fig_hamming_111.png'

pklfile = 'mdf_FN80cysaCA2.pkl'
ot = pd.read_pickle(pklfile)
pklfile = 'mdf_FN28cysaCA2.pkl'
ot2 = pd.read_pickle(pklfile)
pklfile = 'mdf_FN28cysaCA9.pkl'
ot3 = pd.read_pickle(pklfile)

naive28 = ot2[ot2['C28_0']>0].drop(['C28P2CA2','C28P3CA2','C28P5CA2','C28P7CA2'],axis=1)
naive80 = ot[ot['C80_0']>0].drop(['C80P2','C80P3','C80P5','C80P7'], axis=1)

ot = ot.drop(['C80_0'], axis=1)
ot2 = ot2.drop(['C28_0'], axis=1)
ot3 = ot3.drop(['C28_0'], axis=1)

ot = ot.rename(columns={"C80P2": "2", "C80P3": "3", "C80P5": "5", "C80P7": "7"})
ot2 = ot2.rename(columns={"C28P2CA2": "2", "C28P3CA2": "3", "C28P5CA2": "5", "C28P7CA2": "7"})
ot3 = ot3.rename(columns={"C28P2CA9": "2", "C28P3CA9": "3", "C28P5CA9": "5", "C28P7CA9": "7"})

# hamming distance of all variants in winning population
Fn80p2 = ot[ot['2']>0].nlargest(494, '2').drop(['3','5','7'], axis=1)
Fn80p3 = ot[ot['3']>0].nlargest(205, '3').drop(['2','5','7'], axis=1)
Fn80p5 = ot[ot['5']>0].nlargest(169, '5').drop(['2','3','5'], axis=1)
Fn80p7 = ot[ot['7']>0].nlargest(180, '7').drop(['2','3','5'], axis=1)

Fn28p2 = ot2[ot2['2']>0].nlargest(199, '2').drop(['3','5','7'], axis=1)
Fn28p3 = ot2[ot2['3']>0].nlargest(12, '3').drop(['2','5','7'], axis=1)
Fn28p5 = ot2[ot2['5']>0].nlargest(429, '5').drop(['2','3','7'], axis=1)
Fn28p7 = ot2[ot2['7']>0].nlargest(443, '7').drop(['2','3','5'], axis=1)

Fn28p2aCA9 = ot3[ot3['2']>0].nlargest(47, '2').drop(['3','5','7'], axis=1)
Fn28p3aCA9 = ot3[ot3['3']>0].nlargest(83, '3').drop(['2','5','7'], axis=1)
Fn28p5aCA9 = ot3[ot3['5']>0].nlargest(90, '5').drop(['2','3','7'], axis=1)
Fn28p7aCA9 = ot3[ot3['7']>0].nlargest(421, '7').drop(['2','3','5'], axis=1)

df = [naive80, Fn80p2, Fn80p3, Fn80p5, Fn80p7, naive28, Fn28p2, Fn28p3, Fn28p5, Fn28p7, \
     Fn28p2aCA9, Fn28p3aCA9, Fn28p5aCA9, Fn28p7aCA9]

for i in df:
    i['BC_l'] = 0
    i['DE_l'] = 0
    i['FG_l'] = 0
    for j in range(len(i)):
        i['BC_l'].iloc[j] = len(i['BC'].iloc[j])
        i['DE_l'].iloc[j] = len(i['DE'].iloc[j])
        i['FG_l'].iloc[j] = len(i['FG'].iloc[j])
    i = i[(i['BC_l'] < 12) & (i['DE_l'] < 8) & (i['FG_l'] < 10)]

X = 'X'
BC_max, DE_max, FG_max = 9, 6, 8

for i in df:
    i['BC_adj'] = 'A'
    i['DE_adj'] = 'A'
    i['FG_adj'] = 'A'
    i['AA_adj'] = 'A'
    i['AA_adj_arr'] = 'A'
    for j in range(len(i)):
        if i['BC_l'].iloc[j] == BC_max:
            i['BC_adj'].iloc[j] = i['BC'].iloc[j]
        elif i['BC_l'].iloc[j] == BC_max-1:
            i['BC_adj'].iloc[j] = i['BC'].iloc[j][0:4]+X+i['BC'].iloc[j][4:]
        elif i['BC_l'].iloc[j] == BC_max-2:
            i['BC_adj'].iloc[j] = i['BC'].iloc[j][0:3]+2*X+i['BC'].iloc[j][3:]
        if i['DE_l'].iloc[j] == DE_max:
            i['DE_adj'].iloc[j] = i['DE'].iloc[j]   
        elif i['DE_l'].iloc[j] == DE_max-2:
            i['DE_adj'].iloc[j] = i['DE'].iloc[j][0:2]+2*X + i['DE'].iloc[j][2:]
        elif i['DE_l'].iloc[j] == DE_max-3:
            i['DE_adj'].iloc[j] = i['DE'].iloc[j][0:2]+3*X + i['DE'].iloc[j][-1]
        if i['FG_l'].iloc[j] == FG_max:
            i['FG_adj'].iloc[j] = i['FG'].iloc[j]
        elif i['FG_l'].iloc[j] == FG_max-1:
            i['FG_adj'].iloc[j] = i['FG'].iloc[j][0:5]+X+i['FG'].iloc[j][-2:]
        elif i['FG_l'].iloc[j] == FG_max-2:
            i['FG_adj'].iloc[j] = i['FG'].iloc[j][0:5]+2*X+i['FG'].iloc[j][-1:]
    i['AA_adj'] = i['BC_adj']+i['DE_adj']+i['FG_adj']
    for j in range(len(i)):
        i['AA_adj_arr'].iloc[j] = list(i['AA_adj'].iloc[j])
    i['len'] = 0
    for j in range(len(i)):
        i['len'].iloc[j] = len(i['AA_adj'].iloc[j])
        
for i in df:
    i = i[i['len'] == 23 ]
    hamming = None
    hamming = np.empty((len(i),len(i)))
    hamming[:] = np.nan
    for j in range(0,len(i)):
        for k in range(j+1,len(i)):
            hamming[j,k] = sci.hamming(i['AA_adj_arr'].iloc[j], i['AA_adj_arr'].iloc[k]) * 23
    hamming = hamming.flatten()
    hamming = hamming[~np.isnan(hamming)]
    if i.equals(naive80) == True:
        hamN80 = hamming
    elif i.equals(Fn80p2[Fn80p2['len'] == 23]) == True:
        hamFn80p2 = hamming
    elif i.equals(Fn80p3[Fn80p3['len'] == 23]) == True:
        hamFn80p3 = hamming
    elif i.equals(Fn80p5) == True:
        hamFn80p5 = hamming
    elif i.equals(Fn80p7) == True:
        hamFn80p7 = hamming
    elif i.equals(naive28) == True:
        hamN28 = hamming    
    elif i.equals(Fn28p2) == True:
        hamFn28p2 = hamming
    elif i.equals(Fn28p3) == True:
        hamFn28p3 = hamming
    elif i.equals(Fn28p5) == True:
        hamFn28p5 = hamming
    elif i.equals(Fn28p7) == True:
        hamFn28p7 = hamming
    elif i.equals(Fn28p2aCA9) == True:
        hamFn28p2aCA9 = hamming
    elif i.equals(Fn28p3aCA9) == True:
        hamFn28p3aCA9 = hamming
    elif i.equals(Fn28p5aCA9) == True:
        hamFn28p5aCA9 = hamming
    elif i.equals(Fn28p7aCA9[Fn28p7aCA9['len'] == 23]) == True:
        hamFn28p7aCA9 = hamming
    else:
        print('naming is wrong: fix')
    hamming = None

# Generation of best fit line using probability density function
bins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,20,21,22,23]

mu, sigma = norm.fit(np.rint(hamN80))
best_fit_line = norm.pdf(bins, mu, sigma)

mu1, sigma1 = norm.fit(np.rint(hamFn80p2))
best_fit_line1 = norm.pdf(bins, mu1, sigma1)

mu2, sigma2 = norm.fit(np.rint(hamFn80p3))
best_fit_line2 = norm.pdf(bins, mu2, sigma2)

mu3, sigma3 = norm.fit(np.rint(hamFn80p5))
best_fit_line3 = norm.pdf(bins, mu3, sigma3)

mu4, sigma4 = norm.fit(np.rint(hamFn80p7))
best_fit_line4 = norm.pdf(bins, mu4, sigma4)
######
mu5, sigma5 = norm.fit(np.rint(hamN28))
best_fit_line5 = norm.pdf(bins, mu5, sigma5)

mu6, sigma6 = norm.fit(np.rint(hamFn28p2))
best_fit_line6 = norm.pdf(bins, mu6, sigma6)

mu7, sigma7 = norm.fit(np.rint(hamFn28p3))
best_fit_line7 = norm.pdf(bins, mu7, sigma7)

mu8, sigma8 = norm.fit(np.rint(hamFn28p5))
best_fit_line8 = norm.pdf(bins, mu8, sigma8)

mu9, sigma9 = norm.fit(np.rint(hamFn28p7))
best_fit_line9 = norm.pdf(bins, mu9, sigma9)

#####

mu10, sigma10 = norm.fit(np.rint(hamFn28p2aCA9))
best_fit_line10 = norm.pdf(bins, mu10, sigma10)

mu11, sigma11 = norm.fit(np.rint(hamFn28p3aCA9))
best_fit_line11 = norm.pdf(bins, mu11, sigma11)

mu12, sigma12 = norm.fit(np.rint(hamFn28p5aCA9))
best_fit_line12 = norm.pdf(bins, mu12, sigma12)

mu13, sigma13 = norm.fit(np.rint(hamFn28p7aCA9))
best_fit_line13 = norm.pdf(bins, mu13, sigma13)

# Generation of Figure
fig, ax = plt.subplots(figsize=(5,5), dpi=150)

plt.plot(bins, best_fit_line5, color="black", linewidth=2) # naive28
plt.plot(bins, best_fit_line10, color="#1f77b4", linewidth=4) # Fn28p2aCA9
plt.plot(bins, best_fit_line11, color="#ff7f0e", linewidth=4)
plt.plot(bins, best_fit_line12, color="#2ca02c", linewidth=4)
plt.plot(bins, best_fit_line13, color="#d62728", linewidth=4)

plt.plot(bins, best_fit_line, color="black", linewidth=2, linestyle='dashed') # naive80
plt.plot(bins, best_fit_line1, color="#1f77b4", linewidth=2, linestyle='dashed') 
plt.plot(bins, best_fit_line2, color="#ff7f0e", linewidth=2, linestyle='dashed')
plt.plot(bins, best_fit_line3, color="#2ca02c", linewidth=2, linestyle='dashed')
plt.plot(bins, best_fit_line4, color="#d62728", linewidth=2, linestyle='dashed')

plt.plot(bins, best_fit_line6, color="#1f77b4", linewidth=2) # Fn28p2aCA2
plt.plot(bins, best_fit_line7, color="#ff7f0e", linewidth=2)
plt.plot(bins, best_fit_line8, color="#2ca02c", linewidth=2)
plt.plot(bins, best_fit_line9, color="#d62728", linewidth=2)

plt.xlim([0, 23])
plt.ylim([0, 0.2])
plt.xticks(fontsize=14)
plt.yticks([0,0.05,0.1,0.15,0.2],fontsize=14)
plt.xlabel('Hamming Distance', fontsize=16)
plt.ylabel('Normalized Frequency', fontsize=16)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

#plt.legend()
plt.show()
fig.savefig(save, bbox_inches = 'tight')


