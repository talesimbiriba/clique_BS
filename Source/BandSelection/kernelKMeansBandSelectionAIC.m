function [ selectedBandsIdx, Cluster_elem] = kernelKMeansBandSelectionAIC(M,kbw,lambda)
%function [ selectedBandsIdx ] = kernelKMeansBandSelection(M, kbw)
%   M is the endmember matrix
%   kbw is the Gaussian kernel bandwidth



[L,R] = size(M);

if nargin <3
    lambda = 0.5*R;
end

K = computeKernelMatrix(M,kbw);

%maxClusterSize = floor(0.4*L);
%maxClusterSize = floor(0.2*L);
maxClusterSize =45;

count = 1;

aic = 0;
aic_past=1e100;
while (count < maxClusterSize),
    [~, e] = Global_Kernel_K_Means_Clustering(K,ones(L,1),count,'nothing','Fast_Global_KKMeans');
    aic = lambda*count + e;
    if aic>aic_past
        break;
    end
    aic_past = aic;
    count = count+1;
end

Nc = count -1;

[Cluster_elem,Clustering_error] = Global_Kernel_K_Means_Clustering(K,ones(L,1),Nc,'nothing','Fast_Global_KKMeans');


selectedBandsIdx = zeros(Nc,1);
% for over the clusters
for k=1:Nc
    dataIdx = find(Cluster_elem==k);
    Nk = length(dataIdx);
    %initialize zero distance vector
    d = zeros(Nk,1);
    % for over all points in cluster k    
    for n=1:Nk
        % compute the sqr distance of k(.,m_n) and m_k
        d(n) = K(dataIdx(n),dataIdx(n)) -2*sum(K(dataIdx(n),dataIdx))/Nk +  sum(sum(K(dataIdx,dataIdx)))/(Nk*Nk);
    end
    [~,bsIdx]= min(d);
    selectedBandsIdx(k) = dataIdx(bsIdx);
end

selectedBandsIdx = sort(selectedBandsIdx);




end

