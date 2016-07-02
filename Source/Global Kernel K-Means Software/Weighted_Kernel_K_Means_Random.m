function [Cluster_elem,Clustering_error]=Weighted_Kernel_K_Means_Random(K,Dataset_Weights,Clusters,Display,Restarts)
%
%[Cluster_elem,Clustering_error]=Weighted_Kernel_K_Means_Random(K,Dataset_Weights,Total_clusters,Display,Restarts)
%
%This function implements the multiple random restarts Kernel K-Means algorithm as described in  
%G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN. 
%It allows the user to run both the weighted and the non-weighted versions of the Kernel K-Means algorithm.
%
%This function calls the Weighted_Kernel_K_Means function.
%
%Function Inputs
%===============
%
%K is the kernel matrix of the dataset. It must be a positive definite
%square matrix (Gram matrix) in order to guarantee algorithm convergence.
%If K is not positive definite the algorithms may still converge though.
%
%Dataset_Weights is a column vector containing the weight of each datapoint. It 
%is used to run the weighted version of Kernel K-Means. By setting all
%weights equal to 1 the non-weighted version is run. The weights must be positive numbers.
%
%Clusters is the number of clusters.
%
%Display when set to 'nutshell' prints information only about each restart
%and when set to 'details' also prints information about each iteration of Kernel K-Means.
%Any other value results in no printing.
%
%Restarts is the number of random restarts.
%
%Function Outputs
%================
%
%Cluster_elem is a matrix containing the final partitioning of every restart. The clusters are indexed 1,...,Clusters.
%
%Clustering_error is a row vector containing the final clustering error of every restart.
%
%
%Courtesy of G. Tzortzis


%Dataset size.
Data_num=size(K,1);

%Store the partitioning of the dataset for each random restart.
Cluster_elem=zeros(Data_num,Restarts);

%Store the clustering error for each random restart.
Clustering_error=inf(1,Restarts);

for i=1:Restarts
    
%Random initialization of cluster centers.
Init_centers=randperm(Data_num);
Cluster_elem(Init_centers(1:Clusters),i)=1:Clusters;

if strcmp(Display,'details')
    fprintf('\nThe zero clustering error of the first iteration of Kernel K-Means is not valid.\n');
    fprintf('True clustering error is shown after iteration 2. This is due to our implementation.\n');
end

[Cluster_elem(:,i),Clustering_error(i)]=Weighted_Kernel_K_Means(Cluster_elem(:,i),K,Dataset_Weights,Clusters,Display);

if strcmp(Display,'details') || strcmp(Display,'nutshell')
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    fprintf('Restart %d Clusters=%d Clustering Error=%g\n',i,length(unique(Cluster_elem(:,i))),Clustering_error(i));
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
end
end
