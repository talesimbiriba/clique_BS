function [ kbw, lambda ] = searchForBestParametersKKM(y,M,a, lambdas, kbws )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rmseMatrix = zeros(length(kbws),length(lambdas));

[L,N] = size(y);
Nt=floor(0.1*N);
if Nt==0 
    Nt=1;
end
yt = y(:,1:Nt);
at = a(:,1:Nt);



sRMSE = 1e100;
sIdx = [0 0];
i=1;
for kbw=kbws
    j=1;
    for lambda=lambdas
        
        [kkmBS] = kernelKMeansBandSelectionAIC(M, kbw,lambda);
        
        yr=yt(kkmBS,:);
        Mr=M(kkmBS,:);
        a_kkm = tskHype(yr, Mr,[],[],kbw); 
        
        rmseMatrix(i,j) = RMSEAndSTDForMatrix(at,a_kkm);
        if rmseMatrix(i,j) < sRMSE
            sRMSE = rmseMatrix(i,j);
            sIdx = [i j];
        end
        j=j+1;
    end
    i=i+1;
end
kbw = kbws(sIdx(1));
lambda = lambdas(sIdx(2));
end

