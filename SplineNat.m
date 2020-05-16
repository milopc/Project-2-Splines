function [S,b] = SplineNat(vecX, vecY)
%SplineNat This function does the Spline method to interpolate a finite
%amount of points. It is defined in intervals using Hermite polinomials. It
%is a natural Spline, this means an extra condition will be imposed; the
%second derivative of the Spline at the first point and at the last point
%will be set to zero.
%   INPUT:
%vecX: Vector that contains the first coordinates of all the points that
%will be interpolates.
%vecY: Vector that contains the second coordinates of all the points that
%will be interpolated.
%   OUTPUT:
%S: Associated matrix to the system that has to be solved in order to
%obtain the slopes associated to each point that will be interpolated.
%b: Independent coefficients vector associated to the system that has to be
%solved.
%   FUNCION:
n = length(vecX); %This is an auxiliar variable that will help with notation. It is the number of points that will be interpolated.
S = zeros(n,n); %This is the matrix initialized with zeros. For the nature of the matrix that will be worked on, it makes sense that it is a zero matrix.
b = zeros(n,1);%This is the independent coefficients vector initialized with zeros.
%In the next four lines, the S(1,1), S(1,2), S(n-2,n-3), S(n-2,n-2) values
%are assigned manually, it is more efficient than using a loop to define them.
S(1,1)=2;
S(1,2)=1;
S(n,n-1)=1;
S(n,n)=2;
%The vectors h and q contain the absolute change from one coordinate to the
%immediate next on the vecX and vecY respectively.
h = vecX(2:n)-vecX(1:n-1); %Absolute changes in vecX
q = vecY(2:n)-vecY(1:n-1); %Absolute changes in vecY
%This for loop assigns the main diagonal, the entrances to the right and to the left of main diagonal of the S matrix. Also the
%independent coefficients vector.
for i=2:n-1
    S(i,i)=2*(h(i-1)+h(i));
    b(i)=3*((q(i)*h(i-1)/h(i))+(q(i-1)*h(i)/h(i-1)));
    S(i,i-1)=h(i);
    S(i,i+1)=h(i-1);
end
%The first and the last entrances of the coefficient vector have to be
%modified because this process is a Complete Spline.
b(1)=3*q(1)/h(1); 
b(n)=3*q(n-1)/h(n-1);
end

