function [a] = drawAbundacesWithMaxValue(R,numOfPixels, maxAbundance);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = zeros(R,numOfPixels);

for i=1:numOfPixels
    tmp = rand(R,1);
    tmp = tmp/sum(tmp);
    while (max(tmp)>maxAbundance)
        tmp = rand(R,1);
        tmp = tmp/sum(tmp);
    end
    a(:,i) = tmp;
end


end

