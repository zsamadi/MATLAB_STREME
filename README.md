# Matlab STREME

Matlab reimplementation of STREME (Sensitive, Thorough, Rapid, Enriched Motif Elicitation) algorithm which discovers ungapped motifs that are enriched in input sequences. Please refere to algorithm's [official page](https://meme-suite.org/meme/doc/streme.html) for more details. 

The main function for running the algorithm is runStreme. You can run the algorithm on your dataset by runStreme('myFile.fasta'), where myFile.fasta is the name of fasta format file for input sequences. This will run the algorithm with the default settings, listed in the following.

You can change these parameters  by including them in the input parameters, for example,  runStreme('myFile.fasta',W=8, mkvOrder=1, rvp=true) runs the algorithm with motif length set to 8,  Markov order set to 1 and reverse complement set to true. Output results are saved in matStreme.txt file in the output folder.

An example general case of running the algorithm is

```
runStreme('myFile.fasta',NEVAL=25, NREF=4, nRefIter=20, patience=3, evalue=false, nmotifs=0, rvp=true, ...
    mkvOrder=0, W=6,threshold=0.01, hFract=0.1, alphabet='ABCDEFGHIJKLMNO');

```
  
  | input       | Description | 
| :---        |    :----:   |  
| **'exDataFa.fasta'**      | input filename for fasta format data.       | 
| **NEVAL=25**  | number of initial Evaluated seeds.        | 
|  **NREF=4**        |    number of refined seeds.  |
|  **nRefIter= 20**   |    number of iterations for nestled enrichment.  |
|  **patience=3**    |    number consecutive nonsignificant motifs.   |
|  **evalue=false**   |    using evalue instead of pvalue for significance test.  |
|  **nmotifs=0**      |    number of motifs we are looking for, overwrites threshold if it's greater than 0.  |
|  **rvp=true**        |   reverse path.  |
 |  **mkvOrder=0**    |    Markov order for background and negative data generation.  |
|  **W=6**             | length of the motif.   |
|  **threshold=0.01**  |  threshold for motif pvalue obtained from hold-out testing data.   |
|  **hFract=0.1**     |   Fraction of data as hold-out.  |
|  **alphabet='ACGT'**   |     alphabet for input data.  |
