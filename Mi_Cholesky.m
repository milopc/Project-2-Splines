function [L] = Mi_Cholesky(A)
%Mi_Cholesky This function makes the Cholesky factorization process in
%order to express a symmetric and positive-definite matrix A as the product
%of a lower triangular matrix and its traspose. A=L*L'.
%   INPUT:
%A: Symmetric positive-definite matrix
%   OUTPUT:
%L: Lower triangular matrix such that A = L*L'.
n = length(A(:,1));
L = eye(n);
for i = 1:n-1
    L(i,i) = sqrt(A(i,i));
    L(i+1:n,i) = A(i+1:n,i)/L(i,i);
    A(i+1:n,i+1:n) = A(i+1:n,i+1:n)-L(i+1:n,i)*L(i+1:n,i)';
end
L(n,n) = sqrt(A(n,n));
end