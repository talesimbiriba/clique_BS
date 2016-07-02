function [ KM ] = computeKernelMatrix(M, kbw)
%[ KM ] = computeKernelMatrix(M, kbw)
%   Detailed explanation goes here


[L,R] = size(M);

Q=eye(R)/kbw^2;
MQM=M*Q*M';
dMQM = diag(MQM);
KM = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));


KM = (KM + KM')/2;


end

