%% Session 6: 
% Done by Surya Vamshi Penumarthi and Vikas Kurapati
clear all
close all
clc

%% Stage 1: Numerical Integration
x = linspace(0,2,10);
F = x;
I = trapz(F,x);
fprintf('Numerical Integration of x from 0 to 2 is %d.\n',I);

x = 0:0.1:2;
y = 0:0.1:2;
[X,Y] = meshgrid(x,y);
F = cos(X).*sin(Y);
I = trapz(y,trapz(x,F,2));
fprintf('Numerical Integration of cos(x) sin(y) from 0 to 2 is %d.\n',I);

%% Stage 2: Homogeneous Dirichlet Boundary Conditions
Lx = 1;
Ly = 1;

Nx = Lx/2;
Ny = Ly/2;

% 3 source points
x_s = [0.1, 0.5, 0.8];
y_s = [0.1, 0.5, 0.8];

[X_S,Y_S] = meshgrid(x_s,y_s);

lambda = 1.5;
w = 30;

x = linspace(0,1,101); 
y = linspace(0,1,101);
T1 = zeros(length(x));
[X,Y] = meshgrid(x,y);

i = 1;
for x = 0:0.01:1  % Loop over the domain x
    j = 1;
    for y = 0:0.01:1  % Loop over the domain y
        G = 0;
        for m = 1:10 % number of terms considered in Green'S function series expansion
            for n = 1:10
                t1 = 1/((m*pi/Lx)^2 + (n*pi/Ly)^2);
                t2 = sin(m*pi*x/Lx)*sin(m*pi*X_S/Lx)/Nx;
                t3 = sin(n*pi*y/Ly)*sin(n*pi*Y_S/Ly)/Ny;
                G = G + t1*t2*t3;
            end
        end
        T1(i,j) =  w*trapz(y_s,trapz(x_s,G,2))/lambda;
        j = j + 1;
    end
    i = i + 1;
end

figure('Name', 'Surface Plot of the Temperature 1')
surf(X,Y,T1);
title('Temperature Plot with three sources')
subtitle('All boundary conditions are dirichlet')
xlabel('x');
ylabel('y');
zlabel('T');
saveas(gcf,'Stage_2_all_homogeneous_Dirichlet_BC_surface.png')

figure('Name', 'Contour Plot of the Temperature 1')
contourf(X,Y,T1)
title('Temperature Plot with three sources')
subtitle('All boundary conditions are dirichlet')
xlabel('x');
ylabel('y');
zlabel('T');
colorbar
saveas(gcf,'Stage_2_all_homogeneous_Dirichlet_BC_contour.png')

%% Stage 2: Homogeneous Dirichlet (West and South edges) + Homogeneous Neumann (North and East edges) 
Lx = 1;
Ly = 1;

Nx = Lx/2;
Ny = Ly/2;

% 2 source locations
x_s = [0.3,0.7];
y_s = [0.3,0.7];
[X_S,Y_S] = meshgrid(x_s,y_s);
lambda = 1.5;
w = 30;

x = linspace(0,1,101); 
y = linspace(0,1,101);
T2 = zeros(length(x));
[X,Y] = meshgrid(x,y);

i = 1;
for x = 0:0.01:1  
    j = 1;
    for y = 0:0.01:1
        G = 0;
        for m = 1:10
            for n = 1:10
                beta = (2*m-1)*pi/(2*Lx);
                theta = (2*n-1)*pi/(2*Ly);
                t1 = 1/(beta^2 + theta^2);
                t2 = sin(beta*x)*sin(beta*X_S)/Nx;
                t3 = sin(theta*y)*sin(theta*Y_S)/Ny;
                G = G + t1*t2*t3;
            end
        end
        T2(i,j) =  w*trapz(y_s,trapz(x_s,G,2))/lambda;
        j = j + 1;
    end
    i = i + 1;
end

figure('Name', 'Surface Plot of the Temperature 2')
surf(X,Y,T2);
title('Temperature Plot with two sources')
subtitle('West and South BC are dirichlet BC and North and East BC are Neumann')
xlabel('x');
ylabel('y');
zlabel('T');
saveas(gcf,'Stage_2_homogeneous_Dirichlet_Neumann_BC_surface.png')

figure('Name', 'Contour Plot of the Temperature 2')
contourf(X,Y,T2)
title('Temperature Plot with two sources')
subtitle('West and South BC are dirichlet BC and North and East BC are Neumann')
xlabel('x');
ylabel('y');
zlabel('T');
colorbar
saveas(gcf,'Stage_2_homogeneous_Dirichlet_Neumann_BC_contour.png')

%% Stage 3: FD scheme of the domain + source with dirichlet BC on all 4 boundaries 
% Copying the code from previous sessions
Lx = 1;
Ly = 1;
dimX = 101;
dimY = 101;
dx = Lx/dimX;               
dy = Ly/dimY;                
x = linspace(0,Lx,dimX); y = linspace(0,Ly,dimY);
[X,Y] = meshgrid(x,y);
b = zeros(dimX*dimY,1);
index = @(ii,jj) jj + (ii-1)*dimX;
t1 = 1/(dx*dx);
t2 = 1/(dy*dy);
A = zeros(dimX*dimY);

for i = 2 : dimY-1
    for j = 2 : dimX-1
        A(index(i,j),index(i,j)) = -1;
        A(index(i,j),index(i+1,j)) = t1/(2*(t1+t2));
        A(index(i,j),index(i-1,j)) = t1/(2*(t1+t2));
        A(index(i,j),index(i,j+1)) = t2/(2*(t1+t2));
        A(index(i,j),index(i,j-1)) = t2/(2*(t1+t2));
    end
end

%East BC
for i=1:dimY
    A(index(i,dimX),index(i,dimX)) = A(index(i,1),index(i,1)) + 1;
    b(index(i,dimX)) = b(index(i,dimX)) + 4;
end

%West BC
for i=1:dimY
    A(index(i,1),index(i,1)) = A(index(i,1),index(i,1)) + 1;
    b(index(i,1))= b(index(i,1)) + 5;
end

%South BC
for j=1:dimX
    A(index(1,j),index(1,j)) = A(index(1,j),index(1,j)) + 1;
    b(index(1,j)) = b(index(1,j)) + 8;
end

%North BC
for j=1:dimX
    A(index(dimY,j),index(dimY,j)) = A(index(dimY,j),index(dimY,j)) + 1;
    b(index(dimY,j)) = b(index(dimY,j)) + 2;
end

T_finite_diff = A\b;
T_finite_diff = reshape(T_finite_diff,dimY,dimX);

%Adding the solutions of homogeneous Green fuction with three sources and FD scheme with inhomogeneous BC
T_final = T_finite_diff + T1;

figure('Name', 'Surface Plot of the Temperature 3')
surf(X,Y,T_final);
title('Temperature Plot with three sources')
subtitle('FD and GF solution with all boundary conditions are inhomogeneous dirichlet')
xlabel('x');
ylabel('y');
zlabel('T');
saveas(gcf,'Stage_3_FD_GF_Dirichlet_BC_surface.png')

figure('Name', 'Contour Plot of the Temperature 3')
contourf(X,Y,T_final)
title('Temperature Plot with three sources')
subtitle('FD and GF solution with all boundary conditions are inhomogeneous dirichlet')
xlabel('x');
ylabel('y');
zlabel('T');
colorbar
saveas(gcf,'Stage_3_FD_GF_Dirichlet_BC_contour.png')
