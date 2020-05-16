function [L,U] = Mi_LU(A)
%MI_LU Is the function that makes a LU factorization on an square matrix.
%   It returns two matrices, a lower triangular and an upper triangular.
%This solves pont 2 and 3
n = length(A(1,:));%I'm exxtracting the matriz size
L = eye(n);%Im defining a dummy L that I will modify
U = A;%I define a dummy U that has important info and i will modify
for i=1:n-1 %LU factorization process
    L(i+1:n,i) = U(i+1:n,i)/U(i,i);
    U(i+1:n,i+1:n) = U(i+1:n,i+1:n) - (U(i+1:n,i)*U(i,i+1:n))/U(i,i);
    %This proces changes the submatrix on which the algorithm will work in
    %the next iteration
    U(i+1:n,i) = zeros(n-i,1);
    %This step gradually converts the entrances under the main diagonal of
    %the U matrix to ceros
end

