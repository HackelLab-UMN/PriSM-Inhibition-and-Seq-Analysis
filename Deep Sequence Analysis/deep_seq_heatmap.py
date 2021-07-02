import multiprocessing
import numpy as np
from functools import partial
import pandas as pd
pd.options.mode.chained_assignment = None
from sklearn import preprocessing
from joblib import dump, load

# For generation of enrichment heatmap values in Fig 6D. 
# Values were used in Excel to generate heatmap. 
# Python script was run in JupyterLab


pklfile = 'mdf_FN80cysaCA2.pkl'
ot1 = pd.read_pickle(pklfile)
pklfile = 'mdf_FN28cysaCA2.pkl'
ot2 = pd.read_pickle(pklfile)
pklfile = 'mdf_FN28cysaCA9.pkl'
ot3 = pd.read_pickle(pklfile)

# Generate sequences of naive libraries only
naive28 = (ot2[ot2['C28_0']>0].drop(['C28P2CA2','C28P3CA2','C28P5CA2','C28P7CA2'],axis=1)).rename(columns={'C28_0': "naive"})
naive80 = (ot1[ot1['C80_0']>0].drop(['C80P2','C80P3','C80P5','C80P7'], axis=1)).rename(columns={'C80_0': "naive"})
naiveC = naive28.append(naive80)

# Rename columns names to match
ot1 = ot1.rename(columns={"C80_0": "naive","C80P2": "2", "C80P3": "3", "C80P5": "5", "C80P7": "7"})
ot2 = ot2.rename(columns={"C28_0": "naive", "C28P2CA2": "2", "C28P3CA2": "3", "C28P5CA2": "5", "C28P7CA2": "7"})
ot3 = ot3.rename(columns={"C28_0": "naive", "C28P2CA9": "2", "C28P3CA9": "3", "C28P5CA9": "5", "C28P7CA9": "7"})

# Only use sequences that are present in library
ot1['Sum'] = ot1['2']+ot1['3']+ot1['5']+ot1['7']
ot1 = ot1[ot1['Sum'] > 0]
ot2['Sum'] = ot2['2']+ot2['3']+ot2['5']+ot2['7']
ot2 = ot2[ot2['Sum'] > 0]
ot3['Sum'] = ot3['2']+ot3['3']+ot3['5']+ot3['7']
ot3 = ot3[ot3['Sum'] > 0]

# Only use sequences that are present in naive libraries
naive28['Sum'] = naive28['naive']
naive28 = naive28[naive28['Sum']>0]
naive80['Sum'] = naive80['naive']
naive80 = naive80[naive80['Sum']>0]
naiveC['Sum'] = naiveC['naive']
naiveC = naiveC[naiveC['Sum']>0]

df = [ot1, ot2, ot3, naive28, naive80, naiveC]

# determine length of each loop to properly determine frequency
for i in df:
    i['BC_l'] = 0
    i['DE_l'] = 0
    i['FG_l'] = 0
    for j in range(len(i)):
        i['BC_l'].iloc[j] = len(i['BC'].iloc[j])
        i['DE_l'].iloc[j] = len(i['DE'].iloc[j])
        i['FG_l'].iloc[j] = len(i['FG'].iloc[j])
    i = i[(i['BC_l'] < 12) & (i['DE_l'] < 8) & (i['FG_l'] < 10)]

# Insert X for sequences that have less than the max loop length
X = 'X'
BC_max, DE_max, FG_max = 9, 6, 8

for i in df:
    i['BC_adj'] = 'A'
    i['DE_adj'] = 'A'
    i['FG_adj'] = 'A'
    i['AA_adj'] = 'A'
    for j in range(len(i)):
        if i['BC_l'].iloc[j] == BC_max:
            i['BC_adj'].iloc[j] = i['BC'].iloc[j]
        elif i['BC_l'].iloc[j] == BC_max-1:
            i['BC_adj'].iloc[j] = i['BC'].iloc[j][0:2]+X+i['BC'].iloc[j][2:]
        elif i['BC_l'].iloc[j] == BC_max-2:
            i['BC_adj'].iloc[j] = i['BC'].iloc[j][0:2]+2*X+i['BC'].iloc[j][2:]
        if i['DE_l'].iloc[j] == DE_max:
            i['DE_adj'].iloc[j] = i['DE'].iloc[j]   
        elif i['DE_l'].iloc[j] == DE_max-2:
            i['DE_adj'].iloc[j] = i['DE'].iloc[j][0:2]+2*X + i['DE'].iloc[j][2:]
        elif i['DE_l'].iloc[j] == DE_max-3:
            i['DE_adj'].iloc[j] = i['DE'].iloc[j][0:2]+3*X + i['DE'].iloc[j][-1]
        if i['FG_l'].iloc[j] == FG_max:
            i['FG_adj'].iloc[j] = i['FG'].iloc[j]
        elif i['FG_l'].iloc[j] == FG_max-1:
            i['FG_adj'].iloc[j] = i['FG'].iloc[j][0:6]+X+i['FG'].iloc[j][-1:]
        elif i['FG_l'].iloc[j] == FG_max-2:
            i['FG_adj'].iloc[j] = i['FG'].iloc[j][0:5]+2*X+i['FG'].iloc[j][-1:]
    i['AA_adj'] = i['BC_adj']+i['DE_adj']+i['FG_adj']
    
    i['len'] = 0
    for j in range(len(i)):
        i['len'].iloc[j] = len(i['AA_adj'].iloc[j]) 
for i in df:
    i = i[i['len'] == 23 ]
    
    

# One hot encode loop sequences
def encode_seq(encoder,paratope,axis):
    paratope=np.array(list(paratope))
    one_encode=encoder.transform(paratope.reshape(-1,1))  # [ [0,0,1], [0,1,0], ...]   
    return one_encode.flatten() 

for i in df:
    otu_table = i

    AAlist=np.array(list("ACDEFGHIKLMNPQRSTVWXYZ"))         #makes list of 20 amino acids to one hot encode
    encoder=preprocessing.OneHotEncoder(sparse=False)   #sparse=False means will return 2D array instead of sparse matrix
    encoder.fit(AAlist.reshape(-1,1))
    dump(encoder,'one_hot_encoder.joblib')    #saves one hot encoder as joblib
    enc_OH=load('one_hot_encoder.joblib')

    #pass encoder into function that will one-hot encode every row, this helps not have to refit encoder for every sequence
    encode_seq_partial=partial(encode_seq,enc_OH) #or maybe encoder? #enc_OH

    #For every row (axis=1), pass in AAsequence and return one_hot encoded sequence 
    otu_table.loc[:,'One_Hot']=otu_table['AA_adj'].apply(encode_seq_partial,axis=1)

    otu_table.loc[:,'OH_Sum'] = otu_table['One_Hot'] * otu_table['Sum'] #changed rpi to 0 instead of 3

ot1 = ot1[ot1['len'] == 23]
ot2 = ot2[ot2['len'] == 23]
ot3 = ot3[ot3['len'] == 23]
naive28 = naive28[naive28['len'] == 23]
naive80 = naive80[naive80['len'] == 23]

def main(campaign):
    if campaign == 'Fn80':
        ot = ot1
        library = 'Fn80'
    elif campaign == 'Fn28aCA2':
        ot = ot2
        library = 'Fn28aCA2'
    elif campaign == 'Fn28aCA9':
        ot = ot3
        library = 'Fn28aCA9'
    else:
        print('invalid')

    ot_sum = ot['OH_Sum'].sum()
    ot_sum=ot_sum.reshape(23,22)
    ot_sum=ot_sum.transpose()
    ot_freq = ot_sum / ot_sum.sum(axis=0)

    freq1 = ot_freq

    F = freq1[4]
    W = freq1[18]
    Y = freq1[20]
    P = freq1[12]
    M = freq1[10]
    I = freq1[7]
    L = freq1[9]
    V = freq1[17]
    A = freq1[0]
    G = freq1[5]
    C = freq1[1]
    S = freq1[15]
    T = freq1[16]
    N = freq1[11]
    Q = freq1[13]
    D = freq1[2]
    E = freq1[3]
    H = freq1[6]
    K = freq1[8]
    R = freq1[14]
    Z = freq1[21]
    X = freq1[19]
    freq_reordered = np.stack([F,W,Y,P,M,I,L,V,A,G,C,S,T,N,Q,D,E,H,K,R,Z,X])
    freq_reordered
    np.savetxt("./fig6d_" + library + "_freq_reordered.csv", freq_reordered, delimiter=",", fmt="%s")  
    
if __name__ == '__main__':
    main()