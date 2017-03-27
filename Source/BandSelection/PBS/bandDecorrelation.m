function [ bands ] = bandDecorrelation(R, varepsilon )
%function [ bands ] = bandDecorrelation(R, varepsilon )
%  
%  R is the LxN matrix with N pixels with L bands each.
%  varepsilon is the correlation threshold.
%	
% Author: Tales Imbiriba
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
%   Date: December 2016.


[L,n] = size(R);

bands(1) = 1;
bcount = 1;

for i=2:L
    sid = SID(R(i,:), R(bands(bcount),:));
    
    if sid > varepsilon
        bcount = bcount +1;
        bands(bcount) = i;        
    end
    
end



end

