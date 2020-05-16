function [P] = Eval_Hermite(w,a,c)
%Eval_Hermite Evaluates a point on the Hermite polinomial
%   It takes 4 parameters: w is a number, it will be evaluated the Hermite
%   polinomial; a is the first coordinate vector, c is the Hermite
%   polynomial coefficients.
%   It returns the Hermite polinomial evaluated on the point w.
P = c(1) + (w-a(1))*(c(2) + (w-a(1))*(c(3)+c(4)*(w-a(2))));

end

