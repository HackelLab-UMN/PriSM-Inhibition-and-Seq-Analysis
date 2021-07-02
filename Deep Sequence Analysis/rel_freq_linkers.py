import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kde
import seaborn as sns
import matplotlib.ticker as tkr

# Figure 5C
# Relative frequency across different linker length campaigns
# Python script was run in JupyterLab

def main(campaign):
    pklfile = 'mdf_FN80cysaCA2.pkl'
    ot1 = pd.read_pickle(pklfile)
    pklfile = 'mdf_FN28cysaCA2.pkl'
    ot2 = pd.read_pickle(pklfile)
    pklfile = 'mdf_FN28cysaCA9.pkl'
    ot3 = pd.read_pickle(pklfile)

    ot1 = ot1.drop(['C80_0'], axis=1)
    ot2 = ot2.drop(['C28_0'], axis=1)
    ot3 = ot3.drop(['C28_0'], axis=1)

    ot1 = ot1.rename(columns={"C80P2": "2", "C80P3": "3", "C80P5": "5", "C80P7": "7"})
    ot2 = ot2.rename(columns={"C28P2CA2": "2", "C28P3CA2": "3", "C28P5CA2": "5", "C28P7CA2": "7"})
    ot3 = ot3.rename(columns={"C28P2CA9": "2", "C28P3CA9": "3", "C28P5CA9": "5", "C28P7CA9": "7"})

    ot1 = ot1[~((ot1['2']==0) & (ot1['3']==0) & (ot1['5']==0) & (ot1['7']==0))]
    ot2 = ot2[~((ot1['2']==0) & (ot2['3']==0) & (ot2['5']==0) & (ot2['7']==0))]
    ot3 = ot3[~((ot1['2']==0) & (ot3['3']==0) & (ot3['5']==0) & (ot3['7']==0))]

    if campaign == 'Fn80aCA2':
        ot = ot1
        lib = 'Fn80aCA2'
    elif campaign == 'Fn28aCA2':
        ot = ot2
        lib = 'Fn28aCA2'
    elif campaign == 'Fn28aCA9':
        ot = ot3
        lib = 'Fn28aCA9'
    else:
        print('invalid')
        
    arr = ot[['2','3','5','7']].to_numpy()
    freq = arr / arr.sum(axis=0)
    x = np.repeat([[2,3,5,7]],[len(arr)], axis=0)


    # Generation of figure
    fig, ax = plt.subplots(figsize=(4,4), dpi=150)

    for i in range(len(x)):
        plt.plot(x[i],freq[i])
    if lib == 'Fn80aCA2':
        plt.plot(x[17],freq[17],color= '#2ca02c',linewidth=4) # green
        plt.plot(x[12],freq[12],color= '#d62728',linewidth=4) # red
        plt.plot(x[3],freq[3],color= '#ff7f0e',linewidth=4) # orange
        plt.plot(x[4],freq[4],color= "#1f77b4",linewidth=4) # blue
    elif lib == 'Fn28aCA9':
        plt.plot(x[1],freq[1],color= "#ff7f0e",linewidth=4)
        plt.plot(x[2],freq[2],color= '#d62728',linewidth=4)
        plt.plot(x[3],freq[3],color= "#1f77b4",linewidth=4)
        plt.plot(x[4],freq[4],color= '#2ca02c',linewidth=4)
    else:
        plt.plot(x[1],freq[1],color= "#ff7f0e",linewidth=4)
        plt.plot(x[2],freq[2],color= '#d62728',linewidth=4)
        plt.plot(x[3],freq[3],color= "#1f77b4",linewidth=4)
        plt.plot(x[4],freq[4],color= '#2ca02c',linewidth=4)

    plt.xticks([2, 3, 5,7], fontsize=16)
    plt.yticks([0,0.2,0.4], fontsize=16)  #plt.yticks([0,0.5,1])
    plt.ylim([0,0.4])
    ax.set_yticks([0,0.2,0.4]) 
    #ax.set_ylabel('Relative Frequency', labelpad=10)
    ax.margins(x=0, y=0)
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    ax.tick_params(axis='y', direction='out', color='black')
    ax.tick_params(axis='x', length=0)

    plt.tight_layout()

    plt.show()
    #fig.savefig('fig_5C_'+lib+'.png', bbox_inches = "tight")
    # blue orange green red
    # default colors ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
    #    '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
if __name__ == '__main__':
    main()