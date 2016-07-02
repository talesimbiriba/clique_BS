function [Cluster_elem,Clustering_error]=Weighted_Global_Kernel_K_Means_with_CMM(K,Dataset_Weights,Total_clusters,Display,Convits,Top_k,b_frac)
%
%[Cluster_elem,Clustering_error]=Weighted_Global_Kernel_K_Means_with_CMM(K,Dataset_Weights,Total_clusters,Display,Convits,Top_k,b_frac)
%
%This function implements the Global Kernel K-Means with CMM algorithm as described in  
%G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN.
%The CMM algorithm is used to speed up the original Global Kernel K-Means algorithm.
%A Gaussian CMM in feature space is implemented by this code
%which allows the user to run both the weighted and the non-weighted versions of Global Kernel K-Means with CMM.
%
%This function calls the Weighted_Kernel_K_Means and Weighted_Convex_Mixture_Model_Topk functions.
%
%Function Inputs
%===============
%
%K is the kernel matrix of the dataset. It must be a positive definite
%square matrix (Gram matrix) in order to guarantee algorithm convergence.
%If K is not positive definite the algorithms may still converge though.
%
%Dataset_Weights is a column vector containing the weight of each datapoint. It
%is used to run the weighted version of Global Kernel K-Means. By setting all
%weights equal to 1 the non-weighted version is run. The weights must be positive numbers.
%
%Total_clusters is the number of clusters.
%
%Display when set to 'nutshell' prints information only about Global Kernel K-Means
%and when set to 'details' also prints information about each iteration of Kernel K-Means.
%Any other value results in no printing.
%
%Convits is the number of consecutive iterations that the Top_k exemplars in the CMM method
%must remain the same in order the CMM method to converge (usually a value of 1000 suffices).
%
%Top_k is the number of exemplars returned by the CMM method.
%
%b_frac is used to set the b parameter of the CMM method to the value b=b_frac*b0.
%
%Function Outputs
%================
%
%Cluster_elem is a column vector containing the final partitioning of the dataset. The clusters are indexed 1,...,Total_clusters.
%
%Clustering_error is the value of the Kernel K-Means objective function corresponding to the final partitioning.
%
%
%Courtesy of G. Tzortzis

%Dataset size.
Data_num=size(K,1);

%Calculate the empirical distribution of the data.
Data_distribution=Dataset_Weights/sum(Dataset_Weights);

%%%%%%%%%%%%%%FIND GOOD EXEMPLARS TO SPEED UP GLOBAL KK-MEANS%%%%%%%%%%%%%%
if strcmp(Display,'details') || strcmp(Display,'nutshell')
    fprintf('Searching for exemplars\n');
end

%S contains the pairwise Euclidean distances in feature space.
S=-2*K+repmat(diag(K),1,Data_num)+repmat(diag(K)',Data_num,1);

%Calculate the empirical value of the b parameter (b0) for the CMM method.
b=(Data_num*(-log(Data_distribution)'*Data_distribution))/(sum(S,2)'*Data_distribution);
b=b_frac*b;

%Calculate the S matrix for the CMM method.
S=exp(-b*S);

%Cluster the dataset using weighted convex mixture models.
%We use a Gaussian convex mixture model.
[Top_k_idx,Likelihood]=Weighted_Convex_Mixture_Model_Topk(S,Data_distribution',Convits,Top_k);
Top_k_idx=Top_k_idx';
clear S;

if strcmp(Display,'details') || strcmp(Display,'nutshell')
    fprintf('Found exemplars. CMM likelihood=%g\n',Likelihood);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store the optimal error when searching for 1,2,...,Total_clusters.
Best_error=inf(1,Total_clusters);

%Store the assignment of points to clusters corresponding to the optimal solution for 1,2,...,Total_clusters.
Best_clusters=ones(Data_num,Total_clusters);

% Find one cluster solution.
[Best_clusters(:,1),Best_error(1)]=Weighted_Kernel_K_Means(Best_clusters(:,1),K,Dataset_Weights,1,Display);

% Find 2,...,Total_Clusters solutions.
for m=2:Total_clusters
   
    %Solve the m clustering problem by trying all the exemplars found by the CMM method as possible initializations for the m-th cluster.
    for n=Top_k_idx
        Cluster_elem=Best_clusters(:,m-1);
        Cluster_elem(n)=m;
                
       if strcmp(Display,'details') || strcmp(Display,'nutshell')
            fprintf('\n\nSearching for %d clusters.Intitially placing %dth cluster at point %d',m,m,n);
       end
        
       %Find the solution with m clusters.
       [Cluster_elem,Clustering_error]=Weighted_Kernel_K_Means(Cluster_elem,K,Dataset_Weights,m,Display);
       
       if strcmp(Display,'details') || strcmp(Display,'nutshell')
            fprintf('\nFinal Clustering error=%g\n',Clustering_error);
       end
               
        %Keep the best solution with m clusters i.e. the one with lowest clustering error.
        if Best_error(m)>Clustering_error
            Best_error(m)=Clustering_error;
            Best_clusters(:,m)=Cluster_elem;
        end
    end
   
    if size(unique(Best_clusters(:,m)),1)<m
        error('Not able to find more than %d clusters\n',m-1);
    end
end

%Keep as final solution of Global Kernel K-Means with CMM the one with the lowest clustering error (usually its the one with Total_clusters).
[Clustering_error,Clusters]=min(Best_error);
Cluster_elem=Best_clusters(:,Clusters);

if strcmp(Display,'details') || strcmp(Display,'nutshell')
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    fprintf('Best fit:%d clusters with Clustering Error=%g\n',Clusters,Clustering_error);
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
end