#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=24
#SBATCH --mem=200g
#SBATCH -p amdlarge
#SBATCH --job-name="Initial_Prism"
#SBATCH --mail-type=ALL  


module load usearch/11.0_64bit

# merge reads
usearch -fastq_mergepairs "FN28cys0-0"*"_R1_001.fastq" -reverse "FN28cys0-0"*"_R2_001.fastq" -fastqout "FN28cys0_0.fastq" -report "FN28cys0_0.txt" -fastq_minlen 50

usearch -fastq_mergepairs "FN28cys2AAZ0-4-aCA2"*"_R1_001.fastq" -reverse "FN28cys2AAZ0-4-aCA2"*"_R2_001.fastq" -fastqout "FN28cys2PaCA2.fastq" -report "FN28cys2PaCA2.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN28cys3"*"aCA2"*"_R1_001.fastq" -reverse "FN28cys3"*"aCA2"*"_R2_001.fastq" -fastqout "FN28cys3PaCA2.fastq" -report "FN28cys3PaCA2.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN28cys5AAZ0-4-aCA2"*"_R1_001.fastq" -reverse "FN28cys5AAZ0-4-aCA2"*"_R2_001.fastq" -fastqout "FN28cys5PaCA2.fastq" -report "FN28cys5PaCA2.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN28cys7AAZ0-4-aCA2"*"_R1_001.fastq" -reverse "FN28cys7AAZ0-4-aCA2"*"_R2_001.fastq" -fastqout "FN28cys7PaCA2.fastq" -report "FN28cys7PaCA2.txt" -fastq_minlen 50

usearch -fastq_mergepairs "FN28cys2"*"aCA9"*"_R1_001.fastq" -reverse "FN28cys2"*"aCA9"*"_R2_001.fastq" -fastqout "FN28cys2PaCA9.fastq" -report "FN28cys2PaCA9.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN28cys3"*"aCA9"*"_R1_001.fastq" -reverse "FN28cys3"*"aCA9"*"_R2_001.fastq" -fastqout "FN28cys3PaCA9.fastq" -report "FN28cys3PaCA9.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN28cys5"*"aCA9"*"_R1_001.fastq" -reverse "FN28cys5"*"aCA9"*"_R2_001.fastq" -fastqout "FN28cys5PaCA9.fastq" -report "FN28cys5PaCA9.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN28cys7"*"aCA9"*"_R1_001.fastq" -reverse "FN28cys7"*"aCA9"*"_R2_001.fastq" -fastqout "FN28cys7PaCA9.fastq" -report "FN28cys7PaCA9.txt" -fastq_minlen 50

usearch -fastq_mergepairs "FN80cys0-0"*"_R1_001.fastq" -reverse "FN80cys0-0"*"_R2_001.fastq" -fastqout "FN80cys0_0.fastq" -report "FN80cys0_0.txt" -fastq_minlen 50

usearch -fastq_mergepairs "FN80cys2"*"aCA2"*"_R1_001.fastq" -reverse "FN80cys2"*"aCA2"*"_R2_001.fastq" -fastqout "FN80cys2PaCA2.fastq" -report "FN80cys2PaCA2.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN80cys3AAZ0-4-aCA2"*"_R1_001.fastq" -reverse "FN80cys3AAZ0-4-aCA2"*"_R2_001.fastq" -fastqout "FN80cys3PaCA2.fastq" -report "FN80cys3PaCA2.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN80cys5AAZ0-4-aCA2"*"_R1_001.fastq" -reverse "FN80cys5AAZ0-4-aCA2"*"_R2_001.fastq" -fastqout "FN80cys5PaCA2.fastq" -report "FN80cys5PaCA2.txt" -fastq_minlen 50
usearch -fastq_mergepairs "FN80cys7AAZ0-4-aCA2"*"_R1_001.fastq" -reverse "FN80cys7AAZ0-4-aCA2"*"_R2_001.fastq" -fastqout "FN80cys7PaCA2.fastq" -report "FN80cys7PaCA2.txt" -fastq_minlen 50



#Align reads for RPIs
usearch -search_pcr2 "FN28cys0_0.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys0_0_align.fastq"

usearch -search_pcr2 "FN28cys2PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys2PaCA2_align.fastq"
usearch -search_pcr2 "FN28cys3PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys3PaCA2_align.fastq"
usearch -search_pcr2 "FN28cys5PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys5PaCA2_align.fastq"
usearch -search_pcr2 "FN28cys7PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys7PaCA2_align.fastq"

usearch -search_pcr2 "FN28cys2PaCA9.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys2PaCA9_align.fastq"
usearch -search_pcr2 "FN28cys3PaCA9.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys3PaCA9_align.fastq"
usearch -search_pcr2 "FN28cys5PaCA9.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys5PaCA9_align.fastq"
usearch -search_pcr2 "FN28cys7PaCA9.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN28cys7PaCA9_align.fastq"

usearch -search_pcr2 "FN80cys0_0.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN80cys0_0_align.fastq"

usearch -search_pcr2 "FN80cys2PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN80cys2PaCA2_align.fastq"
usearch -search_pcr2 "FN80cys3PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN80cys3PaCA2_align.fastq"
usearch -search_pcr2 "FN80cys5PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN80cys5PaCA2_align.fastq"
usearch -search_pcr2 "FN80cys7PaCA2.fastq" -fwdprimer TCTCTGACTATTTCTTGGGAC -revprimer ATTGATGCTGATTGGGTTTGA -strand plus -fastqout "FN80cys7PaCA2_align.fastq"


#Filter reads, removing any sequences with >1 error or >0 N's
usearch -fastq_filter "FN28cys0_0_align.fastq" -fastaout "FN28cys0_0_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 

usearch -fastq_filter "FN28cys2PaCA2_align.fastq" -fastaout "FN28cys2PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN28cys3PaCA2_align.fastq" -fastaout "FN28cys3PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN28cys5PaCA2_align.fastq" -fastaout "FN28cys5PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN28cys7PaCA2_align.fastq" -fastaout "FN28cys7PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 

usearch -fastq_filter "FN28cys2PaCA9_align.fastq" -fastaout "FN28cys2PaCA9_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN28cys3PaCA9_align.fastq" -fastaout "FN28cys3PaCA9_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN28cys5PaCA9_align.fastq" -fastaout "FN28cys5PaCA9_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN28cys7PaCA9_align.fastq" -fastaout "FN28cys7PaCA9_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 

usearch -fastq_filter "FN80cys0_0_align.fastq" -fastaout "FN80cys0_0_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 

usearch -fastq_filter "FN80cys2PaCA2_align.fastq" -fastaout "FN80cys2PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN80cys3PaCA2_align.fastq" -fastaout "FN80cys3PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN80cys5PaCA2_align.fastq" -fastaout "FN80cys5PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 
usearch -fastq_filter "FN80cys7PaCA2_align.fastq" -fastaout "FN80cys7PaCA2_filter.fasta" -fastq_maxee 1 -fastq_maxns 0 


#find uniques of each sub population
usearch -fastx_uniques FN28cys0_0_filter.fasta -fastaout FN28cys0_0_unique.fasta -sizeout -relabel C28_0-

usearch -fastx_uniques FN28cys2PaCA2_filter.fasta -fastaout FN28cys2PaCA2_unique.fasta -sizeout -relabel C28P2CA2-
usearch -fastx_uniques FN28cys3PaCA2_filter.fasta -fastaout FN28cys3PaCA2_unique.fasta -sizeout -relabel C28P3CA2-
usearch -fastx_uniques FN28cys5PaCA2_filter.fasta -fastaout FN28cys5PaCA2_unique.fasta -sizeout -relabel C28P5CA2-
usearch -fastx_uniques FN28cys7PaCA2_filter.fasta -fastaout FN28cys7PaCA2_unique.fasta -sizeout -relabel C28P7CA2-

usearch -fastx_uniques FN28cys2PaCA9_filter.fasta -fastaout FN28cys2PaCA9_unique.fasta -sizeout -relabel C28P2CA9-
usearch -fastx_uniques FN28cys3PaCA9_filter.fasta -fastaout FN28cys3PaCA9_unique.fasta -sizeout -relabel C28P3CA9-
usearch -fastx_uniques FN28cys5PaCA9_filter.fasta -fastaout FN28cys5PaCA9_unique.fasta -sizeout -relabel C28P5CA9-
usearch -fastx_uniques FN28cys7PaCA9_filter.fasta -fastaout FN28cys7PaCA9_unique.fasta -sizeout -relabel C28P7CA9-

usearch -fastx_uniques FN80cys0_0_filter.fasta -fastaout FN80cys0_0_unique.fasta -sizeout -relabel C80_0-

usearch -fastx_uniques FN80cys2PaCA2_filter.fasta -fastaout FN80cys2PaCA2_unique.fasta -sizeout -relabel C80P2-
usearch -fastx_uniques FN80cys3PaCA2_filter.fasta -fastaout FN80cys3PaCA2_unique.fasta -sizeout -relabel C80P3-
usearch -fastx_uniques FN80cys5PaCA2_filter.fasta -fastaout FN80cys5PaCA2_unique.fasta -sizeout -relabel C80P5-
usearch -fastx_uniques FN80cys7PaCA2_filter.fasta -fastaout FN80cys7PaCA2_unique.fasta -sizeout -relabel C80P7-

#merge filtered reads for same population to calculate ZOTU's
cat FN*cys*PaCA*_filter.fasta > all_filter.fasta

#find uniques of merged populationa
usearch -fastx_uniques all_filter.fasta -fastaout all_unique.fasta -sizeout -relabel Uniq

#denoise unique reads from population, can change minsize for speed
usearch -unoise3 all_unique.fasta -zotus all_zotus.fasta 

#find sequences that match design - create fasta file with sequences that match my design
#>abby1
#ATGNNNATG
#>abby2
#NNNNNKATG
usearch -usearch_global all_zotus.fasta -db designNs.fasta -strand plus -id 1.0 -matched all_zotus_match.fasta

#fill out OTU table for each of the sub populations by the unique file for that sub population
usearch -otutab FN28cys0_0_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys0_0_otutab.txt -threads 10

usearch -otutab FN28cys2PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys2PaCA2_otutab.txt -threads 10
usearch -otutab FN28cys3PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys3PaCA2_otutab.txt -threads 10
usearch -otutab FN28cys5PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys5PaCA2_otutab.txt -threads 10
usearch -otutab FN28cys7PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys7PaCA2_otutab.txt -threads 10

usearch -otutab FN28cys2PaCA9_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys2PaCA9_otutab.txt -threads 10
usearch -otutab FN28cys3PaCA9_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys3PaCA9_otutab.txt -threads 10
usearch -otutab FN28cys5PaCA9_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys5PaCA9_otutab.txt -threads 10
usearch -otutab FN28cys7PaCA9_unique.fasta -zotus all_zotus_match.fasta -otutabout FN28cys7PaCA9_otutab.txt -threads 10

usearch -otutab FN80cys0_0_unique.fasta -zotus all_zotus_match.fasta -otutabout FN80cys0_0_otutab.txt -threads 10

usearch -otutab FN80cys2PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN80cys2PaCA2_otutab.txt -threads 10
usearch -otutab FN80cys3PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN80cys3PaCA2_otutab.txt -threads 10
usearch -otutab FN80cys5PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN80cys5PaCA2_otutab.txt -threads 10
usearch -otutab FN80cys7PaCA2_unique.fasta -zotus all_zotus_match.fasta -otutabout FN80cys7PaCA2_otutab.txt -threads 10


#merge OTU tables across sub-populations
usearch -otutab_merge FN28cys0_0_otutab.txt,FN28cys2PaCA2_otutab.txt,FN28cys3PaCA2_otutab.txt,FN28cys5PaCA2_otutab.txt,FN28cys7PaCA2_otutab.txt -output otutab_FN28cysaCA2.txt
usearch -otutab_merge FN28cys0_0_otutab.txt,FN28cys2PaCA9_otutab.txt,FN28cys3PaCA9_otutab.txt,FN28cys5PaCA9_otutab.txt,FN28cys7PaCA9_otutab.txt -output otutab_FN28cysaCA9.txt
usearch -otutab_merge FN80cys0_0_otutab.txt,FN80cys2PaCA2_otutab.txt,FN80cys3PaCA2_otutab.txt,FN80cys5PaCA2_otutab.txt,FN80cys7PaCA2_otutab.txt -output otutab_FN80cysaCA2.txt



