The Global_Kernel_K_Means_Clustering.m file contains the function than enables the user to
choose between four variants of the Kernel K-Means algorithm and set the appropriate parameters for each variant.

The four variants are Global Kernel K-Means, Fast Global Kernel K-Means, Global Kernel K-Means with Gaussian convex mixture models
and Kernel K-Means with multiple random restarts (both weighted and non-weighted versions can be run) as described in 
G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN.

The software folder contains also the following files:
1)Weighted_Kernel_K_Means.m (Kernel K-Means algorithm)
2)Weighted_Fast_Global_Kernel_K_Means.m (Fast Global Kernel K-Means algorithm)
3)Weighted_Global_Kernel_K_Means.m (Global Kernel K-Means algorithm)
4)Weighted_Global_Kernel_K_Means_with_CMM.m (Global Kernel K-Means with Gaussian convex mixture models algorithm)
5)Weighted_Convex_Mixture_Model_Topk.m (Convex mixture model algorithm)
6)Weighted_Kernel_K_Means_Random.m (Kernel K-Means with multiple random restarts algorithm)
7)3circles_dataset.mat (an example dataset)
8)3circles_kernel_matrix.mat (the 3circles_dataset kernel matrix)


To use this software you must call the Global_Kernel_K_Means_Clustering function and set the appropriate parameters.
The syntax of this function is described below and is also available when typing "help Global_Kernel_K_Means_Clustering" in Matlab.

%[Cluster_elem,Clustering_error]=Global_Kernel_K_Means_Clustering(K,Dataset_Weights,Clusters,Display,Algorithm,varargin)
%
%This functions enables the user to choose between four variants of the Kernel K-Means clustering algorithm.
%These variants are the Global Kernel K-Means, Fast Global Kernel K-Means and Global Kernel K-Means
%with Convex Mixture Models algorithms as proposed in 
%G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN
%as well as the standard Kernel K-Means algorithm with multiple random restarts.
%This software allows the user to run both the weighted and the non-weighted versions of the algorithms described in the above paper.
%Weights are important as when set properly can make the Kernel K-Means objective equivalent to graph cuts.
%
%Function Inputs
%===============
%
%K is the kernel matrix of the dataset. It must be a positive definite
%square matrix (Gram matrix) in order to guarantee algorithm convergence. 
%If K is not positive definite the algorithms may still converge though.
%
%Dataset_Weights is a column vector containing the weight of each datapoint. It
%is used to run the weighted version of the chosen Kernel K-Means variant. By setting all
%weights equal to 1 the non-weighted version is run. The weights must be positive numbers.
%
%Clusters is the number of clusters and must be a positive >=1 integer.
%
%Display when set to 'nutshell' prints synoptic information about the chosen algorithm
%and when set to 'details' also prints information about each iteration of Kernel K-Means.
%Any other value results in no printing.
%
%Algorithm is used to choose the variant to run.
%'Fast_Global_KKMeans' is for the Fast Global Kernel K-Means algorithm,
%'Global_KKMeans' is for the Global Kernel K-Means algorithm,
%'Global_KKMeans_CMM' is for the Global Kernel K-Means with convex mixture models algorithm (a Gaussian CMM is implemented) and
%'KKMeans_Random' is for the Kernel K-Means with multiple random restarts algorithm.
%
%When choosing the 'KKMeans_Random' option an additional parameter is required, the number of restarts.
%The number of restarts must be provided as a scalar after the Algorithm option.
%
%When choosing the 'Global_KKMeans_CMM' option three additional parameters are required, the number of
%iterations (Convits) that the exemplars must remain the same in order the CMM method to converge (usually a value of 1000
%will suffice), the number of returned exemplars (Top_k) by the CMM method and a number (b_frac) than will define the value
%of the b parameter for the CMM method relative to the empirical value b0 i.e. b=b_frac*b0. Setting b_frac=1 results
%in using the empirical b0 value for the method as done in most of the experiments in the above paper.
%The above parameters must be provided as a positive integer, a positive integer and a positive real respectively
%after the Algorithm option.
%
%Function Outputs
%================
%
%Cluster_elem is a column vector containing the final partitioning of the dataset. The clusters are indexed 1,...,Clusters.
%When the 'KKMeans_Random' option is chosen Cluster_elem contains the final partitioning for every restart.
%
%Clustering_error is the value of the Kernel K-Means objective function corresponding to the final partitioning.
%When the 'KKMeans_Random' option is chosen Clustering_error contains the objective function value for every restart.


Usage Example
=============
Suppose we are given the kernel matrix K of a dataset, and we want to partition the dataset into 5 clusters and assign the
same weight to all datapoints.


Cluster the dataset using Fast Global Kernel K-Means
++++++++++++++++++++++++++++++++++++++++++++++++++++
Call Global_Kernel_K_Means_Clustering(K,ones(size(K,1),1),5,'nutshell','Fast_Global_KKMeans')

Cluster the dataset using Global Kernel K-Means
+++++++++++++++++++++++++++++++++++++++++++++++
Call Global_Kernel_K_Means_Clustering(K,ones(size(K,1),1),5,'nutshell','Global_KKMeans')

Cluster the dataset using Global Kernel K-Means with Gaussian CMM
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Call Global_Kernel_K_Means_Clustering(K,ones(size(K,1),1),5,'nutshell','Global_KKMeans_CMM',1000,10,1.2)

1000 defines Convits, 10 is the number of returned exemplars, b=1.2b0 (b0 the empirical b value)

Cluster the dataset using Kernel K-Means with 100 random restarts
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Call Global_Kernel_K_Means_Clustering(K,ones(size(K,1),1),5,'nutshell','KKMeans_Random',100)


Example dataset
===============
Together with the code you can find attached the files 3_circles_dataset.mat and 3_circles_kernel_matrix.mat. By using the kernel matrix
stored in the second file and setting the parameters as described in  
G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN 
you can reproduce the results for the three rings dataset of this paper.
The first file contains the dataset points which make up the three rings.


For any questions regarding the code do not hesitate to contact me at gtzortzi@cs.uoi.gr

G. Tzortzis