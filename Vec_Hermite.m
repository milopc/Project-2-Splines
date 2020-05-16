function [vecOut] = Vec_Hermite(vecIn,a,c)
%Vec_Hermite This function evaluates a function on a vector and returns the
%vector that contains the evaluation of every coordinate.
n = length(vecIn);
vecOut = zeros(n,1);
for i=1:n
    vecOut(i)=Eval_Hermite(vecIn(i),a,c);
end
end

