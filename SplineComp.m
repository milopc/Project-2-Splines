function [S,b] = SplineComp(M1,Mn,vecX,vecY)
%SplineComp This function does the Spline method to interpolate a finite
%amount of points. It is defined in intervals using Hermite polinomials. It
%is a complete Spline, the user will input the slopes at the first and the
%last points of the polinomial.
%   INPUT:
%-M1: Slope at the first point that will be interpolated.
%-M2: Slope at the last point that will be interpolated.
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
S=zeros(n-2,n-2); %This is the matrix initialized with zeros. For the nature of the matrix that will be worked on, it makes sense that it is a zero matrix.
b=zeros(n-2,1); %This is the independent coefficients vector initialized with zeros.
%The vectors h and q contain the absolute change from one coordinate to the
%immediate next on the vecX and vecY respectively.
h = vecX(2:n)-vecX(1:n-1);
q = vecY(2:n)-vecY(1:n-1);
%This for loop assigns the main diagonal of the S matrix and the
%independent coefficients vector.
for i=2:n-1
    S(i-1,i-1)=2*(h(i-1)+h(i));
    b(i-1)=3*(q(i)*h(i-1)/h(i)+q(i-1)*h(i)/h(i-1));
end
%This for loop assigns the entrances above and under the main diagonal on
%the matrix S
for i=2:n-2
    S(i,i-1)=h(i+1);
    S(i-1,i)=h(i-1);
end
%The first and the last entrances of the coefficient vector have to be
%modified because this process is a Complete Spline.
b(1)=b(1)-h(1)*M1; 
b(n-2)=b(n-2)-h(n-2)*Mn;
end

