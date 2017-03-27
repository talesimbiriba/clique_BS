function [ bands ] = bandPrioritizationBS(R, metric, nb, varepsilon)
%[ bands ] = bandPrioritizationBS(R, metric, nb, varepsilon)
%  
%   R is the LxN endmember matrix, where L is the number of bands and N is
%   the number of pixels.
%   
%   metric selects the metric used to prioritize the bands.
%   metric = 0 to select the variance 
%   metric = 1 to select the entropy
%
%   nb is the number of selected bands
%
%   varepsilon is the minimum SID value alowed between close bands
%
%
%   Author: Tales Imbiriba
%   References: 
%
%   [1] C. I. Chang and K. H. Liu, "Progressive Band Selection of Spectral 
%       Unmixing for Hyperspectral Imagery," in IEEE Transactions on 
%       Geoscience and Remote Sensing, vol. 52, no. 4, pp. 2002-2017, April 
%       2014. doi: 10.1109/TGRS.2013.2257604
%
%   [2] C.-I Chang, S. Wang, K. H. Liu, and C. Lin, “Progressive band 
%       dimensionality expansion and reduction for hyperspectral imagery,”
%       IEEE J.  Sel.  Topics  Appl. Earth  Observat. Remote  Sens., vol.4,
%       no. 3, pp. 591–614,  Sep. 2011.
%
%   Date: December 5, 2016.


% varepsilon is the minimum SID value alowed between close bands
if nargin < 4
    varepsilon = 0.001;
end

decorrBandsIDX = bandDecorrelation(R,varepsilon);

% removing correlated bands
R = R(decorrBandsIDX,:);

% compute scores (variance or entropy) for all remaining bands
bp = prioritizeBands( R, metric );

% allocate memory for nb band indexes
bands = zeros(1, nb);

% select the first band
[~,idx] = max(bp);
bp(idx) = -inf;

bands(1) = decorrBandsIDX(idx);
bcount = 1;

% select the rest of the bands.
while bcount <= nb
    % get the index for the largest score
    [~,idx] = max(bp);
    bp(idx) = -inf;
    
    % save the correct band index (see the code line 12 above)
    bands(bcount) = decorrBandsIDX(idx);
    bcount = bcount + 1;
end

% sort the bands before returning the result
bands = sort(bands);

end

