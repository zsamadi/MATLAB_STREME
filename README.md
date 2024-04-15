# MATLAB STREME

MATLAB reimplementation of the STREME (Sensitive, Thorough, Rapid, Enriched Motif Elicitation) algorithm which discovers ungapped motifs that are enriched in input sequences. Please refere to algorithm's [official page](https://meme-suite.org/meme/doc/streme.html) for more details. 

[data](https://github.com/zsamadi/MATLAB-STREME/tree/main/data)

A sample input DNA sequences with fasta format extracted from [RefSeq Transcripts](https://www.ncbi.nlm.nih.gov/genome/guide/human/) dataset.

[run](https://github.com/zsamadi/MATLAB-STREME/tree/main/run)

This folder contains main execution function for the algorithm. Going into this folder and entering streme without any option in the MATLAB command window will run the algorithm with default settings. The default input data is "exData.fasta" in "data" folder. Output results and related figures will be saved in the "output" folder by default. The  settings for the algorithm along with their description and default values are provided in the following table.

  | input       | Description | 
| :---        |    :----:   |  
|'exData.fasta'      | input filename for fasta format data.       | 
|NEVAL=25  | number of initial Evaluated seeds.        | 
|  NREF=4       |    number of refined seeds.  |
|  nRefIter= 20  |    number of iterations for nestled enrichment.  |
|  patience=3    |    number consecutive nonsignificant motifs.   |
|  evalue=false   |    using evalue instead of pvalue for significance test.  |
|  nmotifs=0      |    number of motifs we are looking for, overwrites threshold if it's greater than 0.  |
|  rvc=true        |   reverse complement.  |
 |  mkvOrder=0    |    Markov order for background and negative control data generation.  |
|  W=6            | length of the motif.   |
|  threshold=0.01  |  threshold for motif pvalue obtained from hold-out testing data.   |
| hFract=0.1     |   Fraction of data as hold-out.  |
|  alphabet='ACGT'  |     alphabet for input data.  |

You can use the following command to modify all possible smore options: 

```
streme('exData.fasta',NEVAL=25, NREF=4, nRefIter=20, patience=3, evalue=false, nmotifs=0, rvc=true, ...
    mkvOrder=0, W=6,threshold=0.01, hFract=0.1, alphabet='ACGT');

```
  
[func](https://github.com/zsamadi/MATLAB-STREME/tree/main/func)

Contains related functions for running streme. This folder should be at the same folder as the run folder or its path should be added to the MATLAB path. 

[output](https://github.com/zsamadi/MATLAB-STREME/tree/main/output)

Default folder for saving the output results. The folder currently contains output results of running the algorithm on the default input of "exData.fasta" file. 

## Questions
Please send questions or possible issues with running the code on your data to Zain Samadi (zainsamadi at ucla.edu). 
