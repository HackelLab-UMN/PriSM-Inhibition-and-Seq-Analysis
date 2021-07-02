from Bio import SeqIO
import pandas as pd
import multiprocessing
import numpy as np


# Translates the DNA sequence into AA sequence

# USearch v11 was used to merge and filter deep sequence reads.
# The zotus and otutab files are .txt files generated from USearch analysis.
# Outfile is a pklfile with DNA and AA sequences and sequence counts.

def translate(dna):
    transdic={"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
              "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
              "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
              "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
              "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
              "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
              "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
              "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
              "TAT":"Y","TAC":"Y","TAA":"Z","TAG":"Z",
              "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
              "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
              "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
              "TGT":"C","TGC":"C","TGA":"Z","TGG":"W",
              "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
              "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
              "GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
    AAseq=[]
    if len(dna)%3!=0:
        return "FRAMESHIFT"
    for i in range(0,len(dna),3):
        AAseq.append(transdic[str(dna[i:i+3])])
    AAseq=''.join(AAseq)
    return AAseq

def main(zotus, otutab, outfile):
	#first part loads zotus and translates to AA
    num_cores = multiprocessing.cpu_count()
    pool=multiprocessing.Pool(processes=num_cores) # to run data in parallel, based on number of cpus available
    DNAseq,name=[],[]
    with open(zotus) as indna:        # zotus = 'zotus.fasta'
        parser= SeqIO.parse(indna,'fasta')       # SeqIO.parse can upload a fasta file and read it
        for i in parser:
            DNAseq.append(str(i.seq))
            name.append(i.id)
    (AAseq)=pool.map(translate,DNAseq)

    name_df=pd.DataFrame([name,DNAseq,AAseq])  # pd.DataFrame sets up table where name, DNAseq, AAseq are headers to columns
    name_df=name_df.transpose()              # rows to columns and columns to rows
    name_df.columns=['Zotu_name','DNA','AA']


    #this part merges name_df with the otu-table
    otu_table=otu_table=pd.read_csv(otutab,sep='\t',header=0)    #reading '\t' tab-delimited file 
                                                                                #header=0 keeps headers from 

    name_df['BC'] = name_df['AA'].str.split('YRITY').str[0]
    name_df['BC'].replace('', np.nan, inplace=True)

    name_df['DE'] = name_df['AA'].str.extract('TVPG(.*)ATIS')
    name_df['DE'].replace('', np.nan, inplace=True)

    name_df['FG'] = name_df['AA'].str.split('TITVYAV').str[-1]
    name_df['FG'].replace('', np.nan, inplace=True)

    name_df['AA'].replace('FRAMESHIFT', np.nan, inplace=True)
    name_df=name_df.dropna() #this removes short reads that lack paratope and frameshifts    

    #pd.merge = amazing! Even though zotus.fasta was in random zotu order, still lined up correctly
    merged_otu=pd.merge(name_df,otu_table,how='left',left_on='Zotu_name',right_on='#OTU ID')   #where to merge
    merged_otu=merged_otu.dropna() #this removes short reads that lack paratope and frameshifts    
    merged_otu.to_pickle(outfile)   # 'merged_df.pkl'

if __name__ == '__main__':
    main()
