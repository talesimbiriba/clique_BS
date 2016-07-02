function [ dictionaryIdx ] = buildDictionaryUsingCoherenceFactorKM( K, nu_max )
% [ dictionaryIdx ] = buildDictionaryUsingCoherenceFactor( K1,K2, epsilon0)
%   
%   This function returns the dictionary index dictionaryIdx constructed
%   using a greed coherence based search method. 
%
%   K is the kernel matrix
%   nu_max is the coherence threshold
%
%   Author: Tales Imbiriba
%

if nargin < 2
    nu_max = 0.1;        
end




dictionaryIdx = 1;

L = size(K,1);

if trace(K) == L 
    for l=2:L
        temp = zeros(size(dictionaryIdx));
        for i=1:length(dictionaryIdx)
            temp(i) = K(l,dictionaryIdx(i));
            %temp(i) =  exp(-0.5* (M(l,:) - M(dictionaryIdx(i),:))*(M(l,:) - M(dictionaryIdx(i),:))'/(kbw^2) );
        end
        if max(abs(temp))<=nu_max
            dictionaryIdx = [dictionaryIdx;l];
        end
    end
else

    for l=2:L
        temp = zeros(size(dictionaryIdx));
        for i=1:length(dictionaryIdx)
            temp(i) = K(l,dictionaryIdx(i))/sqrt(K(l,l)*K(dictionaryIdx(i),dictionaryIdx(i)));
            %temp(i) =  exp(-0.5* (M(l,:) - M(dictionaryIdx(i),:))*(M(l,:) - M(dictionaryIdx(i),:))'/(kbw^2) );
        end
        if max(abs(temp))<=nu_max
            dictionaryIdx = [dictionaryIdx;l];
        end
    end

end


end

