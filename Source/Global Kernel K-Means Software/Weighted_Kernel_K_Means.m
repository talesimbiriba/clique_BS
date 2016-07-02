function [Cluster_elem,Clustering_error,Center_dist]=Weighted_Kernel_K_Means(Cluster_elem,K,Dataset_Weights,Clusters,Display)
%
%[Cluster_elem,Clustering_error,Center_dist]=Weighted_Kernel_K_Means(Cluster_elem,K,Dataset_Weights,Clusters,Display)
%
%This function implements the Kernel K-Means algorithm as described in  
%G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN.
%It allows the user to run both the weighted and the non-weighted versions of Kernel K-Means.
%
%This function is called when running Global Kernel K-Means or 
%Fast Global Kernel K-Means or Global Kernel K-Means with CMM or Kernel K-Means with multiple restarts.
%
%Function Inputs
%===============
%
%Cluster_elem is an initial partitioning of the dataset. It is a column vector
%containing the cluster index of each point. The clusters are indexed
%1,2,...,Clusters.
%
%K is the kernel matrix of the dataset. It must be a positive definite
%square matrix (Gram matrix) in order to guarantee algorithm convergence.
%If K is not positive definite the algorithms may still converge though.
%
%Dataset_Weights is a column vector containing the weight of each
%datapoint. It is used to run the weighted version of Kernel K-Means. By setting all
%weights equal to 1 the non-weighted version is run. The weights must be positive numbers.
%
%Clusters is the number of clusters.
%
%Display when set to 'details' prints information after each iteration.
%Any other value results in no printing.
%
%Function Outputs
%================
%
%Cluster_elem is a column vector containing the final partitioning of the dataset. The clusters are indexed 1,...,Clusters.
%
%Clustering_error is the value of the Kernel K-Means objective function corresponding to the final partitioning..
%
%Center_dist is a column vector containing the distance of each point to its cluster center.
%
%
%Courtesy of G. Tzortzis


%Converge if clustering error difference is less than e.
e=0.000001;

%Store the objective function value.
Clustering_error=0;

%Dataset size.
Data_num=size(K,1);

%Store the distance between points and their cluster center.
Center_dist=zeros(Data_num,1);

Iter=1;

while 1
            
    %Keep the clustering error of the previous iteration.
    Old_clustering_error=Clustering_error;
    Clustering_error=0;
       
    Intra=zeros(Clusters,1);
    Cluster_dist=zeros(Data_num,Clusters);
    
    for i=1:Clusters
        
        %Find the dataset points that belong to cluster i and their weights.
        This_elem=find(Cluster_elem==i);
        Dataset_Weights_This_elem=Dataset_Weights(This_elem);
                
        %Calculate intra-cluster pairwise quantity for cluster i.
        Intra(i)=(K(This_elem,This_elem)*Dataset_Weights_This_elem)'*Dataset_Weights_This_elem;
        Intra(i)=Intra(i)/sum(Dataset_Weights_This_elem)^2;
        
        %Calculate point-cluster quantity between all points and cluster i.
        Cluster_dist(:,i)=K(:,This_elem)*Dataset_Weights_This_elem;
        
        %Calculate the distance of all points to the center of cluster i.
        Cluster_dist(:,i)=(-2*Cluster_dist(:,i)/sum(Dataset_Weights_This_elem))+Intra(i)+diag(K);
        
        %Store the distance of cluster's i points to their cluster center.
        Center_dist(This_elem)=Cluster_dist(This_elem,i);
        
        %Add the contribution of cluster i to the clustering error.
        Clustering_error=Clustering_error+Dataset_Weights_This_elem'*Cluster_dist(This_elem,i);
    end
        
    %Update the assignment of points to clusters by placing each point to the closest center.
    [min_dist,Update_cluster_elem]=min(Cluster_dist,[],2);
    
    if strcmp(Display,'details')
        fprintf('\nKernel K-Means Iteration %d\n',Iter);
        fprintf('Clustering error=%g.\n',Clustering_error);
    end
    
    %Check for convergence.
    if Iter>1
        
        if strcmp(Display,'details')
            fprintf('Clustering error reduced by %g\n',Old_clustering_error-Clustering_error);
        end
        
        if abs(Old_clustering_error-Clustering_error) < e
            break;
        end 
    end
          
    Cluster_elem=Update_cluster_elem;
    
    %Drop empty clusters.
    count=0;
    for i=1:Clusters
        if size(find(Update_cluster_elem==i),1)==0
            tmp=find(Update_cluster_elem>i);
            Cluster_elem(tmp)=Cluster_elem(tmp)-1;
            count=count+1;
            warning('Droping empty cluster');
        end
    end 
    
    %Reduce the number of clusters if some were dropped.
    Clusters=Clusters-count;
    
    Iter=Iter+1;
end

return
