function [A, T_init] = construct_matrix(dimX, dimY, boundary)
% Initialization
dx = 1/(dimX-1);                     % Spatial step in x direction
dy = 1/(dimY-1);                    % Spatial step in y direction

% Define Temperature T vector (Intial values with zero for example)
T_init = zeros(dimX*dimY, 1);

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
%% Adding Boundary conditions
for j = 1:dimY
    A(index(1,j), index(1,j)) = 1;
    T_init(index(1, j)) = boundary.T_D;
end

for i = 2:dimX-1
    if strcmp(boundary.south, 'Dirichlet')
        A(index(i,1), index(i,1)) = 1;
        T_init(index(i,1)) = boundary.T_D;
    else
        A(index(i,1), index(i,1)) = 1.5*boundary.lambda/dy;
        A(index(i,1),index(i,2)) = -2*boundary.lambda/dy;
        A(index(i,1), index(i,3)) = 0.5*boundary.lambda/dy;
        T_init(index(i, 1)) = T_init(index(i, 1)) + boundary.heat_flux/boundary.lambda;
    end
end

for i = 2:dimX-1
    if strcmp(boundary.north, 'Dirichlet')
        A(index(i, dimY), index(i, dimY)) = 1;
        T_init(index(i, dimY)) = boundary.T_D;
    else
        A(index(i,dimY), index(i,dimY)) = boundary.alpha + 1.5*boundary.lambda/dy;
        T_init(index(i, dimY)) = boundary.alpha*boundary.T_inf;
        A(index(i, dimY), index(i, dimY-1)) = -2*boundary.lambda/dy;
        A(index(i, dimY), index(i, dimY-2)) = 0.5*boundary.lambda/dy;
    end
end

for j = 1:dimY
    if strcmp(boundary.east, 'Dirichlet')
        A(index(dimX, j), index(dimX, j)) = 1;
        T_init(index(dimX, j)) = boundary.T_D;
    else
        A(index(dimX, j), index(dimX, j)) = boundary.alpha + 1.5*boundary.lambda/dx;
        T_init(index(dimX, j)) = boundary.alpha*boundary.T_inf;
        A(index(dimX, j), index(dimX-1, j)) = -2*boundary.lambda/dx;
        A(index(dimX, j), index(dimX-2, j)) = 0.5*boundary.lambda/dx;
    end
end

%%-----------------------------------------------------------------------
%% Defining the Source

T_init(index(dimX/2, dimY/2)) = T_init(index(dimX/2, dimY/2)) + 100; %midde of the domain. This defines extra heat flux.
end