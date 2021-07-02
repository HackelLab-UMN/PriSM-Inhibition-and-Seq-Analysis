import numpy as np
import pandas as pd


# Figure 6C
# Determination of relative frequency for each loop length of each campaign
# Script was run in JuypterLab
# Final plot was created in Excel

# Campaigns are defined as otFn80, otFn28aCA2, & otFn28aCA2
def main(c):
    pklfile = 'mdf_FN80cysaCA2.pkl'
    ot = pd.read_pickle(pklfile)
    pklfile = 'mdf_FN28cysaCA2.pkl'
    ot2 = pd.read_pickle(pklfile)
    pklfile = 'mdf_FN28cysaCA9.pkl'
    ot3 = pd.read_pickle(pklfile)

    naive28 = ot2[ot2['C28_0']>0].drop(['C28P2CA2','C28P3CA2','C28P5CA2','C28P7CA2'],axis=1)
    naive80 = ot[ot['C80_0']>0].drop(['C80P2','C80P3','C80P5','C80P7'], axis=1)

    ot = ot.rename(columns={"C80_0": "naive","C80P2": "2", "C80P3": "3", "C80P5": "5", "C80P7": "7"})
    ot2 = ot2.rename(columns={"C28_0": "naive", "C28P2CA2": "2", "C28P3CA2": "3", "C28P5CA2": "5", "C28P7CA2": "7"})
    ot3 = ot3.rename(columns={"C28_0": "naive", "C28P2CA9": "2", "C28P3CA9": "3", "C28P5CA9": "5", "C28P7CA9": "7"})

    otFn80 = ot
    otFn28aCA2 = ot2
    otFn28aCA9 = ot3

    list_df = [otFn80, otFn28aCA2, otFn28aCA9]

    for df in list_df:
        df['BC_l'] = 0
        df['DE_l'] = 0
        df['FG_l'] = 0

    for df in list_df:
        for i in range(len(df)):
            df['BC_l'].iloc[i] = len(df['BC'].iloc[i])
            df['DE_l'].iloc[i] = len(df['DE'].iloc[i])
            df['FG_l'].iloc[i] = len(df['FG'].iloc[i])
        df = df[(df['BC_l'] < 12) & (df['DE_l'] < 8) & (df['FG_l'] < 10)]   

    if c == 'otFn80':
        a = list_df[0]
    elif c == 'otFn28aCA2':
        a = list_df[1]
    elif c == 'otFn28aCA9':
        a = list_df[2]
    else:
        print('invalid')

    df = a
    df = df[df['naive']>0]
    naive = [ [0, 0, len(df[df['FG_l'] == 6])], \
             [len(df[df['BC_l'] == 7]), len(df[df['DE_l'] == 3]), len(df[df['FG_l'] == 7])], \
             [len(df[df['BC_l'] == 8]), len(df[df['DE_l'] == 4]), len(df[df['FG_l'] == 8])], \
             [len(df[df['BC_l'] == 9]), len(df[df['DE_l'] == 6]), 0] ]
    df = a
    df = df[df['2']>0]
    peg2 = [ [0, 0, len(df[df['FG_l'] == 6])], \
            [len(df[df['BC_l'] == 7]), len(df[df['DE_l'] == 3]), len(df[df['FG_l'] == 7])], \
            [len(df[df['BC_l'] == 8]), len(df[df['DE_l'] == 4]), len(df[df['FG_l'] == 8])], \
            [len(df[df['BC_l'] == 9]), len(df[df['DE_l'] == 6]), 0] ] 
    df = a
    df = df[df['3']>0]
    peg3 = [ [0, 0, len(df[df['FG_l'] == 6])], \
            [len(df[df['BC_l'] == 7]), len(df[df['DE_l'] == 3]), len(df[df['FG_l'] == 7])], \
            [len(df[df['BC_l'] == 8]), len(df[df['DE_l'] == 4]), len(df[df['FG_l'] == 8])], \
            [len(df[df['BC_l'] == 9]), len(df[df['DE_l'] == 6]), 0] ]
    df = a
    df = df[df['5']>0]
    peg5 = [ [0, 0, len(df[df['FG_l'] == 6])], \
            [len(df[df['BC_l'] == 7]), len(df[df['DE_l'] == 3]), len(df[df['FG_l'] == 7])], \
            [len(df[df['BC_l'] == 8]), len(df[df['DE_l'] == 4]), len(df[df['FG_l'] == 8])], \
            [len(df[df['BC_l'] == 9]), len(df[df['DE_l'] == 6]), 0] ] 
    df = a
    df = df[df['7']>0]
    peg7 = [ [0, 0, len(df[df['FG_l'] == 6])], \
            [len(df[df['BC_l'] == 7]), len(df[df['DE_l'] == 3]), len(df[df['FG_l'] == 7])], \
            [len(df[df['BC_l'] == 8]), len(df[df['DE_l'] == 4]), len(df[df['FG_l'] == 8])], \
            [len(df[df['BC_l'] == 9]), len(df[df['DE_l'] == 6]), 0] ] 
    
    fn = np.array(naive)/np.sum(naive, axis=0)
    f2 = np.array(peg2)/np.sum(peg2, axis=0)
    f3 = np.array(peg3)/np.sum(peg3, axis=0)
    f5 = np.array(peg5)/np.sum(peg5, axis=0)
    f7 = np.array(peg7)/np.sum(peg7, axis=0)
    
    f_list = [fn, f2, f3, f5, f7]
    f_list = np.array(f_list).reshape(20,3)
    np.savetxt("./loop_length_" + c + ".csv", f_list, delimiter=",", fmt="%s") 
    print('Complete')
    
    
if __name__ == '__main__':
    main()