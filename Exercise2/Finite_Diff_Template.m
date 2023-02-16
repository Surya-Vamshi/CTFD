%% PR - Computational Thermo Fluid Dynamics
%  TUM - Ass. Professorship for Thermo Fluid Dynamics
%  SS 022


clear all; 
close all; 
clc;

%% Numerical code to solve the 2D steady heat equation by Finite Differences
%  Second order accurate in Space

% Geometrical Parameters

l = 1;                      % Length of the domain in x 
h = 1;                      % Length of the domain in y
dimX = 100;                   % Number of nodes in x
dimY = 100;                 % Number of nodes in y 
%nDof =                   % Number of degrees of freedom

dx = 1/(dimX-1);                     % Spatial step in x direction
dy = 1/(dimY-1);                    % Spatial step in y direction

% Defining the point heat source
pointsource = 'off';
            % 1) on
            % 2) off

% if on set
%xheatsource =    % x position
%yheatsource =    % y position
%omega       =  % intensity


% Defining Boundary conditions
% Type
boundary.south = 'Neumann';
boundary.north = 'Robin';
boundary.east  = 'Robin';
boundary.west  = 'Dirichlet';

% Value for  BC


% Thermal conductivity Coefficient

%heat_conduc = 'non_homogenous';
            % 1) homgenous      2) non_homogenous (region with different K)
            % 3) random         4) linear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------
% Define Temperature T vector (Intial values with zero for example)
T_init = zeros(dimX*dimY, 1);

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
% Define the Heat conductivity coefficient values

Lambdamatrix = ones(dimX*dimY, 1);


% switch heat_conduc
%     
%     case 'homogenous'
%         
%     case 'non_homogenous'        
%         
%     case 'random'
%         
%     case 'linear'
% 
% end


%%------------------------------------------------------------------------
%% Constructing computational matrix;
index = @(ii, jj) ii + (jj-1)*dimX;
dx2 = dx*dx;
dy2 = dy*dy;

A = zeros(dimX*dimY);

for i = 2:dimX-1
    for j = 2:dimY - 1
        A(index(i,j), index(i,j)) = -2/(dx2) - 2/dy2;
        A(index(i,j), index(i+1, j)) = 1/(dx2);
        A(index(i,j), index(i-1, j)) = 1/(dx2);
        A(index(i,j), index(i, j+1)) = 1/(dy2);
        A(index(i,j), index(i, j-1)) = 1/(dy2);
    end
end

%%-----------------------------------------------------------------------
%% Defining Boundary conditions
T_inf = 400;
T_D = 600;
alpha = 3;
lambda = 1.5;
heat_flux = 500;

for j = 1:dimY
    A(index(1,j), index(1,j)) = 1;
    T_init(index(1, j)) = T_D;
end

for i = 2:dimX-1
    A(index(i,1), index(i,1)) = 1.5*lambda/dy;
    A(index(i,1),index(i,2)) = -2*lambda/dy;
    A(index(i,1), index(i,3)) = 0.5*lambda/dy;
    T_init(index(i, 1)) = T_init(index(i, 1)) + heat_flux/lambda;
end

for i = 2:dimX-1
    A(index(i,dimY), index(i,dimY)) = alpha + 1.5*lambda/dy;
    T_init(index(i, dimY)) = alpha*T_inf;
    A(index(i, dimY), index(i, dimY-1)) = -2*lambda/dy;
    A(index(i, dimY), index(i, dimY-2)) = 0.5*lambda/dy;
end

for j = 1:dimY
    A(index(dimX, j), index(dimX, j)) = alpha + 1.5*lambda/dx;
    T_init(index(dimX, j)) = alpha*T_inf;
    A(index(dimX, j), index(dimX-1, j)) = -2*lambda/dx;
    A(index(dimX, j), index(dimX-2, j)) = 0.5*lambda/dx;
end

%%-----------------------------------------------------------------------
%% Defining the Source

T_init(index(dimX/2, dimY/2)) = T_init(index(dimX/2, dimY/2)) + 1000; %midde of the domain. This defines extra heat flux.

T_init(index(1, dimY/2)) = T_init(index(1, dimY/2)) + 50; %source at dirichilet boundary. But this fixes temperature

%%------------------------------------------------------------------------      
%% Solving the linear system (use '\' operator)
T_final = A\T_init;

%------------------------------------------------------------------------
% Ploting Results

% convert solution vector into matrix
T_final = reshape(T_final, dimX, dimY);

% Do a surface plot (use 'surf()' function)
surf(X, Y, T_final');
savefig("3D_temperature_plot.fig")

% Do a contour plot (use 'contour()' and 'colorbar' functions)
figure(2)
contourf(X, Y, T_final')
colorbar
savefig("contours.fig")