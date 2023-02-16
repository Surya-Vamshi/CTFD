function T = solveFVM(M, X, Y, boundary, TD, alpha, Tinf, lambda, q_dot_sym)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it
%
% input
% M         Spatial Matrix M
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
%
% output
% T         Temperature field

% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) * size(M, 1);
[dimY, dimX] = size(M);

% B is the right-hand side of the linear system
B = zeros(dimX*dimY, 1);

% Set boundary conditions
if strcmp(boundary.north, 'Dirichlet')
    B(index(1,1:dimX)) = TD.north;
end
if strcmp(boundary.south, 'Dirichlet')
    B(index(dimY,1:dimX)) = TD.south;
end
if strcmp(boundary.east, 'Dirichlet')
    B(index(1:dimY,dimX)) = TD.east;
end
if strcmp(boundary.west, 'Dirichlet')
    B(index(1:dimY,1)) = TD.west;
end


% set up the system matrix A
A = zeros(dimY*dimX);

verbose = 0;
for i = 1:size(M, 1)
    for j = 1:size(M, 2)
        % Fill the system matrix and the right-hand side for node (i,j)
        [A(index(i,j), :), B(index(i,j))] =  stamp(i, j, X, Y, B(index(i,j)), alpha, TD, Tinf, lambda, q_dot_sym, boundary, verbose);
    end
end

% solve the linear system
T(:) = A\B;




