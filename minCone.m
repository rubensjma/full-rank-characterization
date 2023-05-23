function [H,b] = minCone(R)

% R: (2^{n-1} x n) array of points

% V Points defining the set conv(V ).
% R Rays defining the set cone(R).

n = size(R,2);
V = zeros(1,n);
P = Polyhedron('V',V,'R',R);
Haux = P.H;
H = Haux(:,1:end-1);
b = Haux(:,end);