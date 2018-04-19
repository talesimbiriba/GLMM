function [ RMSE, STD] = RMSEAndSTDForMatrix(A, B)
%[ RMSE, STD] = RMSEAndSTDForMatrix(A, B)
%   Compute RMSE and STD for matrices A and B.

[R,N] = size(A);

E = A - B;
RMSE = sqrt(trace(E'*E)/(N*R));
STD = sqrt(var((diag(E'*E))/(R)));

end

