## Usage

From the root directoy, just type make. The executable will be located in bin folder with name runner. 

The application requires the CUDA compiler nvcc with device architecture of 3.5 capability. Thrust library is also required. 

The .cu files are in the /src folder and all the header files are located in /include directory. When compiling a build directory will be automatically created and the exectubale is located in /bin. 

To execute the program:
 
```
#./bin/runner <TFFile> <PatientDataFile> <nTFs> <nGenes> <nSamples> <nBootstraps> <pValue>  > outputFile
./bin/runner data/tfs data/input_expmat_512_clean.txt 1813 20531 512 1 0.0000001 > data/result_3_col_512.txt
```

The data files are not included in the git repo, but is available through download. 

PatienDataFile: https://www.dropbox.com/s/42wg6imiu3efl30/brca-expmat-512-clean.csv?dl=0

TFFile: https://www.dropbox.com/s/b36h2nf1j9fhmgy/tfs.txt?dl=0

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

Because the queue data structure is in shared memory, a careful calculation of its upper bound is essential for parallel execution. The theoretical upper bound is (nSamples - 3) items in the queue. 

## Data processing inequality

Here it expects an input of matrix with size of nTFs * nTFs * nGenes.

nTFs * nTFs * nGenes threads will be launched to prune  the graph. The result is saved in a boolean matrix mask and later applied using nTFs * nGenes threads. 

Another note: interaction with itself is not taken care of, aka. the degenerated triangle is left without special consideration to reduce thread divergence. 

## Bootstrapping

Random sample the same number of columns from the original data, ship it to device.

2 additional matrix initialized on device, one for the sum of all MIs for each bootstrap, and the other for counts of occurence this edge among all bootstraps.

Overall edges and occurence were calculated on device: using parallel reduce with thrust library. 

The 2 bootstrap matrix along with the total edges and occurence are shipped back to host

Poisson model is generated with mean edge calculated from total occurence / total edge, and for each edge, the p-value was calculated based on the occurence of this edge

For each edge, if possion pvalue is less then 0.05(bonferroni corrected), then tf, gene, meanMI, and poisson pvalue were be write out. 

Output: 3 column text file 

## Current performance

For 200 samples, 20531 genes, and 1813 transcription factors. The time taken (including IO) is 64 sec. 

Comparing to the multithreadded Java version running on CPU which takes around 192 sec (4 core) and 649 sec (1 core) to execute 200 samples. 

