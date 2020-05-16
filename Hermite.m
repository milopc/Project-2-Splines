function [c] = Hermite(a,b,m)
%Hermite Calculates the Hermite polinomial coefficients
%   It takes three vectors. a is the first coordinates vector, b is the
%   second coordinates vector and m is the slopes vector associated to the
%   (a1,b1),(a2,b2) points.
c = zeros(4,1);
h = a(2)-a(1);
q = (b(2)-b(1))/h;
c(1) = b(1);
c(2) = m(1);
c(3) = (q-m(1))/h;
c(4) = (m(1) + m(2) - 2*q)/(h^2);
end

