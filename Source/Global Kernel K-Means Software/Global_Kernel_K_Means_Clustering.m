function [Cluster_elem,Clustering_error]=Global_Kernel_K_Means_Clustering(K,Dataset_Weights,Clusters,Display,Algorithm,varargin)
%
%[Cluster_elem,Clustering_error]=Global_Kernel_K_Means_Clustering(K,Dataset_Weights,Clusters,Display,Algorithm,varargin)
%
%This functions enables the user to choose between four variants of the Kernel K-Means clustering algorithm.
%These variants are the Global Kernel K-Means, Fast Global Kernel K-Means and 
%Global Kernel K-Means with Convex Mixture Models algorithms as proposed in 
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
%Clusters is the number of clusters and must be a postive >=1 integer.
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
%
%
%Courtesy of G. Tzortzis

% tic
if nargin<5
    error('Too few input arguments.Type help Global_Kernel_K_Means_Clustering for proper use');
elseif nargin>8 
    error ('Too many input arguments.Type help Global_Kernel_K_Means_Clustering for proper use');
end

%Check if the input arguments have the correct properties.

%Check the kernel matrix.
if length(size(K))~=2
    error('K must be a 2D matrix');
else
    [K_rows,K_columns]=size(K);
    if K_rows~=K_columns
        error('K must be square');
    end
end

%Check the vector with the dataset weights.
if length(size(Dataset_Weights))>2
    error('Dataset_Weights must be a vector');
else
    [W_rows,W_columns]=size(Dataset_Weights);
    if W_rows~=K_rows || W_columns~=1
        error ('Dataset_Weights must be a column vector with as many rows as matrix K');
    elseif ~isempty(find(Dataset_Weights<=0,1))
        error('All weights must be positive numbers in Dataset_Weights');
    end
end

%Check the Clusters variable.
if size(Clusters,1)~=1 || size(Clusters,2)~=1 || length(size(Clusters))>2
    error('Clusters must be a scalar');
elseif Clusters<1 
    error('The number of clusters must be greater or equal to 1');
else
    Clusters=floor(Clusters);
    if K_rows<Clusters
        error('K must have more rows than the number of clusters');
    end
end

%Check Display.
if ~ischar(Display)
    error('Display must be a character string. ''nutshell'' for synoptic printing ''details'' for analytical printing or any other value for no printing');
end

%Check the Algorithm variable.
if strcmp(Algorithm,'Fast_Global_KKMeans') && nargin==5
    [Cluster_elem,Clustering_error]=Weighted_Fast_Global_Kernel_K_Means(K,Dataset_Weights,Clusters,Display);
elseif strcmp(Algorithm,'Global_KKMeans') && nargin==5
    [Cluster_elem,Clustering_error]=Weighted_Global_Kernel_K_Means(K,Dataset_Weights,Clusters,Display);
elseif strcmp(Algorithm,'KKMeans_Random') && nargin==6 && varargin{1}>=1
    Restarts=floor(varargin{1});
    [Cluster_elem,Clustering_error]=Weighted_Kernel_K_Means_Random(K,Dataset_Weights,Clusters,Display,Restarts);
elseif strcmp(Algorithm,'Global_KKMeans_CMM') && nargin==8 && varargin{1}>=1 && varargin{2}>=1 && varargin{3}>0
    Convits=floor(varargin{1});
    Top_k=floor(varargin{2});
    b_frac=varargin{3};
    [Cluster_elem,Clustering_error]=Weighted_Global_Kernel_K_Means_with_CMM(K,Dataset_Weights,Clusters,Display,Convits,Top_k,b_frac);
else
    error('Algorithm must be ''Fast_Global_KKMeans'' or ''Global_KKMeans'' or ''KKMeans_Random'' or ''Global_KKMeans_CMM''. For the ''KKMeans_Random'' choice an additional argument is required, the number of Restarts>=1 while for the ''Global_KKMeans_CMM'' choice three additional arguments (type "help Global_Kernel_K_Means_Clustering" for details)' );
end
% toc;


