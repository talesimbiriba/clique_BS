function [ sid ] = SID( b1, b2 )
%function [ sid ] = SID( b1, b2 )
%   SID - Spectral Information Divergence
%
%   based in the Kullback-Leiber information measure
%   b1 - unormalized band 1
%   b2 - unormalized band 2
%
%   Author: Tales Imbiriba
%   December 2016

p1 = b1./(sum(b1));
p2 = b2./(sum(b2));


sid = sum(p1.*log(p1./p2)) +  sum(p2.*log(p2./p1)) ;


end

