import multiprocessing
import numpy as np
from functools import partial
import pandas as pd
pd.options.mode.chained_assignment = None
from sklearn import preprocessing
from joblib import dump, load

# Used for one hot encoding amino acid sequence to a series of 0s and 1s
# To use, input pklfile and newfile names.

def encode_seq(encoder,paratope,axis):
    paratope=np.array(list(paratope))
    one_encode=encoder.transform(paratope.reshape(-1,1))  # [ [0,0,1], [0,1,0], ...]   
    return one_encode.flatten()   # removes individual lists from array [0,0,1,0,1,0,...] "flattens" whole thing

# sort (ex: 'bind23BC')
# pklfile (ex: 'mdf_bind23BC.pkl')
# newfile (ex: 'oh_mdf_bind23BC.pkl')

def main(pklfile, newfile):
    #load dataframe
    otu_table=pd.read_pickle(pklfile)     #reads merged.pkl file  #'./merged_df.pkl'

    #Create one-hot-encoder for each of the 21 amino acids (including Z) in particular order!!!!
    AAlist = np.array(list("FWYPMILVAGCSTNQDEHKRZX"))         #makes list of 21 amino acids to one hot encode
    encoder = preprocessing.OneHotEncoder(sparse=False)   #sparse = False means will return 2D array instead of sparse matrix
    encoder.fit(AAlist.reshape(-1,1))                 #makes into one row instead of 20 rows x 1
    #joblib is more efficient on objects that carry large numpy arrays internally - can only pickle to disk and not to string
    dump(encoder,'one_hot_encoder.joblib')    #saves one hot encoder as joblib
    enc_OH = load('one_hot_encoder.joblib')
    
    #pass encoder into function that will one-hot encode every row, this helps not have to refit encoder for every sequence
    encode_seq_partial = partial(encode_seq,enc_OH)

    #For every row (axis = 1), pass in AAsequence and return one_hot encoded sequence 
    BC = otu_table['BC']
    DE = otu_table['DE']
    FG = otu_table['FG']
    otu_table.loc[:,'OH_BC'] = otu_table['BC'].apply(encode_seq_partial,axis=1)
    otu_table.loc[:,'OH_DE'] = otu_table['DE'].apply(encode_seq_partial,axis=1)
    otu_table.loc[:,'OH_FG'] = otu_table['FG'].apply(encode_seq_partial,axis=1)

     
    #save otu table
    otu_table.to_pickle(newfile)
    print('Complete')


    
    
    
    

if __name__ == '__main__':
    main()