INFOH304-BIOINF PROJECT: PROTEIN SEQUENCE ALIGNMENT USING THE SMITH-WATERMAN ALGORITHM

To produce a protein sequence alignment using the Smith-Waterman algorithm and/or an exact sequence match, 
all the files are included in the folder 'swNew'.
The database used for the test is 'newE.fasta', containing all the proteins which names start with the letter E. 
To execute the code for a test protein 'P00533.fasta': 
Add in the same folder the database's binary files with the extensions '.pin', '.psq' and '.phr'.
Run the command:  

        make
        
And then:
  
        ./projet ./newE.fasta ./P00533.fasta

The default parameters are:  BLOSUM62, gap_opening=11, gap_extension=1 
To change the parameters run:

        ./projet ./newE.fasta ./P00533.fasta -o IntOpenPenalty -e IntExtensionPenalty -b Blosum
                            



