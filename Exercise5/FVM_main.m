clear all
close all
clc

% SESSION_05

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady and unsteady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables

InitFVM

%% initialize spatial Matrix T

M = zeros(dimY,dimX);

%% set up the mesh

[X, Y] = setUpMesh(M, l, formfunction);

%% Fill matrix A and vector B. Solve the linear system.
solveFVM(M, X, Y, boundary, TD, alpha, Tinf, lambda, q_dot_sym, problem, t_end, theta, solver, no_iter, tol);

