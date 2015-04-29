## Data input

Input is via pure file stream read. 1.3 sec to load a 20000 * 256 matrix. 

Time for loading all gene labels is not measurable. 

Time for loading all TF names is not measurable.

Creating a mapping between the index of TFs and the gene indexes used unordered map on CPU, which might be slow. Could use GPU but need to deal with string loading on GPU memory. 

## Rank all rows

Ranking is done with nRows threads on GPU. Each thread calls the thrust library function with scheduling policy. 

Ranking 20000 * 2560 matrix takes 4.3 sec.

## Constructing null model

Random permuted data were fed into the GPU kernel (this indeed incures data shipping). 100000 pairs of random permuted samples were calculated for mutual information, which gives 100000 mutual information values. They are then sorted on GPU. The outliers (100) were thrown away. The cumulative mass function was computed and the tail is fit with linear regression on CPU. 

## Calculating mutual information 

Kernel with thread dimension of nGenes * nGenes * nSamples is launched. Each pair of mutual information is computed with a BLOCK of threads. Threads in the block are used for counting. A queue data structure is kept in shared memory to record the state of adaptive partitioning.  

4 additional values are kept in the shared memory for counting purpose. Counting is carried out using hard-coded bit operation. Also repeated overwrite of the same data is allowed to reduce thread divergence. 

## Data processing inequality

Here it expects an input of matrix with size of nTFs * nTFs * nGenes.

nTFs * nTFs * nGenes threads will be launched to prune  the graph. The result is saved in a boolean matrix mask and later applied using nTFs * nGenes threads. 

One more thing to think of: if during DPI, two edges are equally small, which one should we cut?

Another note: interaction with itself is not taken care of, aka. the degenerated triangle is left without special consideration to reduce thread divergence. 

So far, the whole process from null model construction and data processing inequality takes around 3 minutes (including around 600 MB disk IO) with sample size of 200, gene number of 20000. Note that the Java multithread version running on 32 core computing cluster takes around 9 min (including disk IO). 

## Bootstrapping

random sample the same number of columns from the original data, ship it to device

2 additional matrix on device, one for the sum of all MIs for each bootstrap, and the other for counts of occurence this edge among all bootstraps

overall edges and occurence were calculated on device

the 2 bootstrap matrix along with the total edges and occurence are shipped back to host

poisson model is generated with mean edge calculated from total occurence / total edge, and for each edge, the p-value was calculated based on the occurence of this edge

for each edge, if possion pvalue is less then 0.05(bonferroni corrected), then tf, gene, meanMI, and poisson pvalue were be write out. 

Output: 4 column text file 