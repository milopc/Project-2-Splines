function [] = Graph_Hermite(Mata,Matb,Matm, colGros)
%Graph_Hermite This function graphs the Hermite Polinomial obtained from
%joining length(Mata(1,:)) intervals.
%   Mata is the matrix whose column i are the the a_j vector of first
%   coordinates of the function at the begining and at the end of the
%   inteval i.
%   Matb is the matrix whose column i are the the b_j vector of second
%   coordinates of the function at the begining and at the end of the
%   inteval i.
%   Matm is the matrix whose columns are the the c_j vector of the slope
%   of the function at the begining and at the end of the inteval i.
%   colGros is the string parameter that will indicate the type of line and
%   the color.
n = length(Mata(1,:));
for i=1:n
    c = Hermite(Mata(:,i),Matb(:,i),Matm(:,i));
    x = linspace(Mata(1,i),Mata(2,i),20);
    y = Vec_Hermite(x,Mata(:,i),c);
    plot(x,y,colGros, 'LineWidth',1.5);
    title ("Resultados gráficos de electrocardiograma")
    xlabel ("Tiempo en milisegundos")
    ylabel ("Voltios")
    hold on
end
end

