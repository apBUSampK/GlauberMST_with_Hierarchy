# GlauberMST with Hierarchical Clustering


## Idea

The idea of the project is to apply more advanced clustering algorithms to GlauberMST in order to resolve the discrepancy of the result of AAMCC and experimental results for light nuclei fragmentation (see [article](https://doi.org/10.3390/particles5010004)), while using hierarchical clustering methods to optimize searching for an optimal clusterization.

Note: this is a prototype heavily based on GlauberMST.


## Implemented

* Algorithm for building c-linkage dendrogram
* Silhouette comparison for finding the best clusterization near the given critical distance (see [Silhouette (clustering)](https://en.wikipedia.org/wiki/Silhouette_(clustering)))

## To-do list

* Analyse the behavior of silhouette-based clustering
* Implement an entropy algorithm for finding the best clusterization (based on number of protons and neutrons in a cluster, see [Entropy (information theory)](https://en.wikipedia.org/wiki/Entropy_(information_theory)))
* Combine entropy and silhouette algorithms
* Try other decision rules
* Integrate into [AAMCC](https://github.com/alexsvetlichnyy/AAMCC)
