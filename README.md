## Data input

Input is via pure file stream read. 1.3 sec to load a 20000 * 256 matrix. 
Time for loading all gene labels is not measurable. 
Time for loading all TF names is not measurable.
Creating a mapping between the index of TFs and the gene indexes used unordered map on CPU, which might be slow. Could use GPU but need to deal with string loading on GPU memory. 

## Rank all rows

Ranking is done with nRows threads on GPU. Each thread calls the thrust library function with scheduling policy. 
Ranking 20000 * 2560 matrix takes 4.3 sec.

## Constructing null model

## Calculating mutual information 

## Data processing inequality

Here it expects an input of matrix with size of nTFs * nTFs * nGenes.

nTFs * nTFs * nGenes threads will be launched to prune  the graph. The result is saved in a boolean matrix mask and later applied using nTFs * nGenes threads. 

## Bootstrapping

