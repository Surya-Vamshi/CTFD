clear all
close all
clc

% SESSION_03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables

InitFVM

%% initialize spatial Matrix T

M = zeros(dimY,dimX);

%% set up the mesh

[X, Y] = setUpMesh(M, l, formfunction);

%% Fill matrix A and vector B. Solve the linear system.

T = solveFVM(M, X, Y, boundary, TD, alpha, Tinf, lambda, q_dot_sym);

%% Make some plots
T = reshape(T, dimY, dimX);

T = [T; flip(T)];
X = [X; flip(X)];
Y = [Y; flip(-Y)];

figure(1)
fig = contourf(X, Y, T);
colormap("jet");
colorbar;
saveas(gcf, "cooling_fin2d.fig")

figure(2)
surf(X, Y, T);
colormap("jet");
colorbar;
saveas(gcf, "cooling_fin3d.fig")