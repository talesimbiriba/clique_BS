function [ RMSE, STD] = RMSEAndSTDForMatrix_bigData(A, B)
%[ RMSE, STD] = RMSEAndSTDForMatrix(A, B)
%   Compute RMSE and STD for matrices A and B.

[R,N] = size(A);

E = A - B;

e_norm=zeros(N,1);
for i=1:N,
    e_norm(i) = E(:,i)'*E(:,i);
end

RMSE = sqrt(sum(e_norm)/(N*R));
STD = sqrt(var(e_norm/R));


% RMSE = sqrt(trace(E'*E)/(N*R));
% STD = sqrt(var((diag(E'*E))/(R)));

end

