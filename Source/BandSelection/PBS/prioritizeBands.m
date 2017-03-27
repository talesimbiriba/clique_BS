function [ bands_p ] = prioritizeBands( R, metric )
%function [ bandsp ] = prioritizeBands( R, metric )
%
%   R - Image matrix Lxn, where L is the number of bands and n is the
%   number of pixels
%
%   metric = 
%              0 - variance
%              1 - entropy
%


[L,n] =  size(R);

if metric == 0
    bands_p = var(R');
elseif metric ==1
%     Rn = R/(sum(sum(R)));
%     bands_p = -sum( (Rn.*log(Rn))' );
    bands_p = zeros(1,L);
    for i=1:L
        bands_p(i) = entropy(R(i,:));
    end
end

end

