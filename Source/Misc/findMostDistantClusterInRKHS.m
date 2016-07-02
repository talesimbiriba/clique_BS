function [ clusterNumber ] = findMostDistantClusterInRKHS( Cluster_elem, K )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


L = size(K,1);
numberOfClusters = max(Cluster_elem);
sqrDists = zeros(numberOfClusters,1);

for k=1:numberOfClusters
    clusterIdx = find(Cluster_elem==k); 
    Nk = length(clusterIdx);
    
    sqrDists(k) = (1/(Nk^2)) * sum(sum(K(clusterIdx,clusterIdx))) -(2/(L*Nk))*sum(sum(K(clusterIdx,:))) + (1/L^2)*sum(sum(K));
    
end

[~,clusterNumber] = max(sqrDists);

end

