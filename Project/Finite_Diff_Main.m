%% Computational Thermo Fluid Dynamics
%  Surya Vamshi Penumarthi and Vikas Kurapati
% Project


clear all;
close all;
clc;

%% Numerical code to solve the 2D steady heat equation by Finite Differences
%  Second order accurate in Space

% Geometrical Parameters

l = 1;                      % Length of the domain in x
h = 1;                      % Length of the domain in y
dimX = 16;                   % Number of nodes in x
dimY = 16;                 % Number of nodes in y
%nDof =                   % Number of degrees of freedom

dx = 1/(dimX-1);                     % Spatial step in x direction
dy = 1/(dimY-1);                    % Spatial step in y direction

%Choose between solver types
% 1) direct
% 2) jacobi
% 3) gauss-siedel
% 4) SOR
% 5) steepest_descent
% 6) conjugate_gradient
% 7) bi_conjugate_gradient
% 8) conjugate_gradient_squared
% 9) gmres
% 10) multigrid_2level
% 11) multigrid_vcycle
% 12) multigrid_wcycle
% 13) multigrid_fcycle

solver = 'multigrid_fcycle';
no_iter = 25000;
tol = 1e-5;

max_level = 3;  % Maximum number of levels for multigrid


% Defining Boundary conditions
% Type
boundary.south = 'Robin';
boundary.north = 'Neumann';
boundary.east  = 'Robin';
boundary.west  = 'Dirichlet';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------
% Create mesh (use 'meshgrid()' function).

x = linspace(0,l,dimX);
y = linspace(0,h, dimY);
[X,Y] = meshgrid(x,y);

% coordinate system should look like this:
%        Y^(0,dimY)
%         |
%         |
%         |
%         |
%         |
%         |
%         |____________________>X
%    (0,0)                   (0,dimX)


%----------------------------------------------------------------------
%% Defining Boundary conditions
boundary.T_inf = 400;
boundary.T_D = 100;
boundary.alpha = 3;
boundary.lambda = 1.5;
boundary.heat_flux = 500;

%% Getting A and T_init
[A, T_init] = construct_matrix(dimX, dimY, boundary);
T_guess = boundary.T_D*ones(size(T_init));

%%------------------------------------------------------------------------
%% Solving the linear system (using various solvers)
switch solver
    case 'direct'
        tic
        T_final = A\T_init;
        time_taken = toc;
        status = 1;

    case 'jacobi'
        tic
        [T_final, status] = solver_jacobi(A, T_init, T_guess, no_iter, tol);
        time_taken = toc;

    case 'gauss-siedel'
        tic
        [T_final, status] = solver_gauss_seidel(A, T_init, T_guess, no_iter, tol);
        time_taken = toc;

    case 'SOR'
        tic
        [T_final, status] = solver_sor(A, T_init, T_guess, no_iter, tol,1.7);
        time_taken = toc;

    case 'steepest_descent'
        tic
        [T_final, status] = solver_steepest_descent(A, T_init, T_guess, no_iter, tol);
        time_taken = toc;

    case 'conjugate_gradient'
        tic
        [T_final, status] = solver_conjugate_gradient(A, T_init, T_guess, no_iter, tol);
        time_taken = toc;

    case 'bi_conjugate_gradient'
        tic
        [T_final, status] = solver_bi_conjugate_gradient(A, T_init, T_guess, no_iter, tol);
        time_taken = toc;

    case 'conjugate_gradient_squared'
        tic
        [T_final, status] = solver_conjugate_gradient_squared(A, T_init, T_guess, no_iter, tol);
        time_taken = toc;

    case 'gmres'
        tic
        [T_final, status] = solver_gmres(A, T_init, T_guess, no_iter, tol, 100);
        time_taken = toc;

    case 'multigrid_2level'
        tic
        [T_final, status] = multigrid_2level(A, T_init, T_guess, dimX, dimY, boundary, no_iter, tol);
        time_taken = toc;

    case 'multigrid_vcycle'
        tic
        [T_final, status] = multigrid_vcycle(A, T_init, T_guess, dimX, dimY, boundary, no_iter, tol, max_level);
        time_taken = toc;

    case 'multigrid_wcycle'
        tic
        [T_final, status] = multigrid_wcycle(A, T_init, T_guess, dimX, dimY, boundary, no_iter, tol, max_level);
        time_taken = toc;

    case 'multigrid_fcycle'
        tic
        [T_final, status] = multigrid_fcycle(A, T_init, T_guess, dimX, dimY, boundary, no_iter, tol, max_level);
        time_taken = toc;

    otherwise
        disp('Please select a solver type')
end

%% Ploting Results if solution converges
if status == 1
    disp(append('Solved steady state problem using ',solver, ' method'))
    fprintf('Time taken to solve = %f seconds\n',time_taken)

    T = reshape(T_final, dimY, dimX)';

%     figure(1)
%     contourf(X, Y, T);
%     title(append('Contour Plot of the Temperature, solved using ',solver, ' method'))
%     colormap("jet");
%     colorbar;
% 
%     figure(2)
%     surf(X, Y, T);
%     title(append('Surface Plot of the Temperature, solved using ',solver, ' method'))
%     colormap("jet");
%     colorbar;
else
    disp(append('Unable to solve steady state problem using ',solver, ' method'))
end
