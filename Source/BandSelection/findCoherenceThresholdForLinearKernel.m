function [ mu_0, bands_Idx ] = findCoherenceThresholdForLinearKernel(M, numOfBands, bs_alg )
%UNTITLED2 Summary of this function goes here
%   Searchs for the mu_0 that produces a dictionary with numOfBands bands.
%   M is the LxR endmember matrix
%   numOfBands is the number of desired bands
%   bs_alg sets the bs algorithm. 
%   use bs_alg = 0 (default) to perform a greedy band selection
%   use bs_alg = 1 to perform a maximum clique band selection
%   

if nargin < 3
    bs_alg = 0;
end


mu_0 = 0.1;
KM = M*M';


% search for the mu parameter starting in mu_0 = 0.1 and incrementing

if bs_alg==0
    nb = size(buildDictionaryUsingCoherenceFactorKM( KM, mu_0 ),1);

    step = 0.1;

    while (nb ~= numOfBands) 

        if (nb< numOfBands) %&& (mu_0+step)<1,
            mu_0 = mu_0+step;
            nb = size(buildDictionaryUsingCoherenceFactorKM( KM, mu_0 ),1);
            while nb>numOfBands
                mu_0 = mu_0 -step;
                step = step/2;
                mu_0 = mu_0 + step; 
                nb = size(buildDictionaryUsingCoherenceFactorKM( KM, mu_0 ),1);
            end
        end

    end

    bands_Idx = buildDictionaryUsingCoherenceFactorKM( KM, mu_0 );

else
    
    % normalize kernel matrix (unity norm matrix)
    KM = KM./sqrt(diag(KM)*diag(KM)');
    
    nb = size(clique_coherence_bandselection( KM, mu_0, [], 1 ),1);

    step = 0.1;

    while (nb ~= numOfBands) 

        if (nb< numOfBands) %&& (mu_0+step)<1,
            mu_0 = mu_0+step;
            if mu_0==1 
                nb = size(KM,1);
            else
                nb = size(clique_coherence_bandselection( KM, mu_0, [], 1 ),1);
            end
            while nb>numOfBands
                mu_0 = mu_0 -step;
                step = step/2;
                mu_0 = mu_0 + step; 
                nb = size(clique_coherence_bandselection( KM, mu_0, [], 1 ),1);
            end
        end

    end

    bands_Idx = clique_coherence_bandselection( KM, mu_0, [], 1 );

    
end



end

