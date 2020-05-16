%       PROJECT  2     COMPUTATIONAL  MATHEMATICS
%           Emiliano Padilla Cardona - 166136
%           Salvador Rodríguez Carrasco - 171943 

%FIRST. Two vectors are defined, the vecX contains the value at which the
%60 equally spaced intervals between 1 and 1000 start. The vecY contains
%the observations at each of the vecX entrances.
vecX = linspace(1,1000,60);
vecY = xlsread("DatosECG.xlsx");


%       PART I. Composite Spline interpolation.
%   I.1: The matrix associated with the Composite Spline interpolation problem
%is obtained. The coefficients vector associated with the problem is
%obtained.
Slope1 = 50; %The imposed slope at the first point.
Slopen = -50; %The imposed slope at the last point.
[S,b]=SplineComp(Slope1,Slopen,vecX,vecY);
%   I.2: 
%The S matrix is diagonally strictly dominant matrix. Also because of how
%the data is provided, it is symmetric. Then it can be concluded that it is
%symmetric and positive-definite. Therefore, it has cholesky factorization.
%The vector that contains the slopes of the polinomial from the second to
%the fifty-ninth point is obtained.
L1 = Mi_Cholesky(S);
y=Solve_LT(L1,b);
M1=Solve_UT(L1',y);
%   I.3: Three matrices will be created
%The first matrix contains on its first row the vecX values except for the 
%last and on its second row the vecX values except for the first. Its 
%columns represent each interval that will be worked with.
Mata=[vecX(1:59);vecX(2:60)];
%The second matrix contains on its first row the vecY values except for the 
%last and on its second row the vecY values except for the first. Its 
%columns represent the observation at the begining and at the end of each
%interval.
Matb=[vecY(1:59)';vecY(2:60)'];
%The third matrix contains on its first row the slopes of the points except
%for the last and on its second row the slopes of the points except for the
%first.
Matm=[Slope1 M1';M1' Slopen];

%   I.4: The graph of the Spline is generated through a function.
Graph_Hermite(Mata, Matb, Matm,'b')




%       PART II. Natural Spline Interpolation.

%   II.1: The matrix associated with the Composite Spline interpolation problem
%is obtained. The coefficients vector associated with the problem is
%obtained.
[T,d]=SplineNat(vecX,vecY);

%   II.2: 
%The slopes vector will be obtained by solving the system defined by the T
%matrix and the d vector with the LU factorization. The T matrix is not
%symmetric so it can't be factorized with Cholesky.
[L2,U2] = Mi_LU(T);
t=Solve_LT(L2,d);
M2=Solve_UT(U2,t);

%   II.3: Three matrices will be created
%The first matrix contains on its first row the vecX values except for the 
%last and on its second row the vecX values except for the first. Its 
%columns represent each interval that will be worked with.
MataN=[vecX(1:59);vecX(2:60)];
%The second matrix contains on its first row the vecY values except for the 
%last and on its second row the vecY values except for the first. Its 
%columns represent the observation at the begining and at the end of each
%interval.
MatbN=[vecY(1:59)';vecY(2:60)'];
%The third matrix contains on its first row the slopes of the points except
%for the last and on its second row the slopes of the points except for the
%first.
MatmN=[M2(1:59)';M2(2:60)'];

%   II.4: The graph of the Spline is generated through a function.
Graph_Hermite(MataN, MatbN, MatmN, '--g')