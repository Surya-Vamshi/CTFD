function    solveFVM(M, X, Y, boundary, TD, alpha, Tinf, lambda, q_dot_sym, problem, t_end, theta, solver, no_iter, tol)


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
% problem   Whether problem is steady or unsteady case
% t_end     Simulation end time for unsteady case
% theta     theta value for the weighted average or θ-method
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
A= sparse(A);

verbose = 0;
for i = 1:size(M, 1)
    for j = 1:size(M, 2)
        % Fill the system matrix and the right-hand side for node (i,j)
        [A(index(i,j), :), B(index(i,j))] =  stamp(i, j, X, Y, B(index(i,j)), alpha, TD, Tinf, lambda, q_dot_sym, boundary, verbose);
    end
end


%% Solving problem
switch problem
    case "steady"
        %% solve the linear system
        % Initial guess for iterative solvers
        T = zeros(dimX*dimY,1);
        status = 0;
        switch solver
            case 'direct'
                tic
                T(:)= A\B;
                time_taken = toc;
                fprintf('Time taken to solve = %f \n',time_taken)
                whos A;
                status = 1;

            case 'jacobi'
                [T(:), status] = solver_jacobi(A,B,T,no_iter,tol);

            case 'gauss-siedel'
                [T(:), status] = solver_gauss_seidel(A,B,T,no_iter,tol);

            case 'SOR'
                [T(:), status] = solver_sor(A,B,T,no_iter,tol,1.7);

            otherwise
                disp('Please select a solver type in InitFVM.m')
        end

        %% Make some plots if solution converges
        if status == 1
            disp(append('Solved steady state problem using ',solver, ' method'))
            
            T = reshape(T, dimY, dimX);

            T = [T; flip(T)];
            X = [X; flip(X)];
            Y = [Y; flip(-Y)];

            figure(1)
            contourf(X, Y, T);
            title(append('Contour Plot of the Temperature, solved using ',solver, ' method'))
            colormap("jet");
            colorbar;

            figure(2)
            surf(X, Y, T);
            title(append('Surface Plot of the Temperature, solved using ',solver, ' method'))
            colormap("jet");
            colorbar;
        else
            disp(append('Unable to solve steady state problem using ',solver, ' method'))
        end

    case "unsteady-explicit"
        % Initial guess of the temperature was taken as Tinf
        T = ones(size(M,1)*size(M,2),1)*Tinf;

        % Applying Dirichlet boundary condition to the initial guess
        if strcmp(boundary.north, 'Dirichlet')
            T(index(1,1:dimX)) = TD.north;
        end
        if strcmp(boundary.south, 'Dirichlet')
            T(index(dimY,1:dimX)) = TD.south;
        end
        if strcmp(boundary.east, 'Dirichlet')
            T(index(1:dimY,dimX)) = TD.east;
        end
        if strcmp(boundary.west, 'Dirichlet')
            T(index(1:dimY,1)) = TD.west;
        end

        % Calculating dt based on stability condition
        s_factor = 0.25; % safety_factor

        x_min = abs(min(X(1, 2:end) - X(1, 1:end-1)));
        y_min = abs(min(Y(2:end, dimX) - Y(1:end-1, dimX)));

        l_c= sqrt(x_min*y_min);

        dtc = s_factor * (l_c^2)/(4*lambda);

        X = [X; flip(X)];
        Y = [Y; flip(-Y)];

        progress = waitbar(0,'Code is runnning...');
        delete cooling_fin2d_unsteady_explicit.gif
        delete cooling_fin3d_unsteady_explicit.gif

        %% Solving the usteady problem
        for t=0:dtc:t_end
            waitbar(t/t_end, progress)

            % Using explicit euler method
            T = T + dtc*(A*T-B) ;

            %% Make some plots
            T_plot = reshape (T, dimY, dimX);
            T_plot  =[T_plot; flip(T_plot)];

            f1=figure(1);
            contourf(X, Y, T_plot);
            title(['Contour Plot of the Temperature at time=',num2str(t),'(sec)'])
            colormap("jet");
            colorbar;
            % exportgraphics(gcf,"cooling_fin2d_unsteady_explicit.gif",'Append',true);
            movegui(f1,'east')

            f2=figure(2);
            surf(X, Y, T_plot);
            title(['Surface Plot of the Temperature at time=',num2str(t),'(sec)'])
            colormap("jet");
            colorbar;
            % exportgraphics(gcf,"cooling_fin3d_unsteady_explicit.gif",'Append',true);
            movegui(f2,'west')
        end

        close(progress)
        disp('Solved unsteady state problem using Explicit Euler method')

    case "unsteady-implicit"

        % Initial guess of the temperature was taken as Tinf
        T = ones(size(M,1)*size(M,2),1)*Tinf;

        % Applying Dirichlet boundary condition to the initial guess
        if strcmp(boundary.north, 'Dirichlet')
            T(index(1,1:dimX)) = TD.north;
        end
        if strcmp(boundary.south, 'Dirichlet')
            T(index(dimY,1:dimX)) = TD.south;
        end
        if strcmp(boundary.east, 'Dirichlet')
            T(index(1:dimY,dimX)) = TD.east;
        end
        if strcmp(boundary.west, 'Dirichlet')
            T(index(1:dimY,1)) = TD.west;
        end

        % Setting dt to a fixed value
        dt = 0.1;

        X = [X; flip(X)];
        Y = [Y; flip(-Y)];

        progress = waitbar(0,'Code is runnning...');
        delete cooling_fin2d_unsteady_implicit.gif
        delete cooling_fin3d_unsteady_implicit.gif

        %% Solving the usteady problem
        for t=0:dt:t_end
            waitbar(t/t_end, progress)

            % Using implicit theta method
            T = (eye(size(A))-theta*dt*A)\(T+(1-theta)*dt*A*T - dt* B);

            %% Make some plots
            T_plot = reshape (T, dimY, dimX);
            T_plot  =[T_plot; flip(T_plot)];

            f1=figure(1);
            contourf(X, Y, T_plot);
            title(['Contour Plot of the Temperature at time=',num2str(t),'(sec)'])
            colormap("jet");
            colorbar;
            % exportgraphics(gcf,"cooling_fin2d_unsteady_implicit.gif",'Append',true);
            movegui(f1,'east')

            f2=figure(2);
            surf(X, Y, T_plot);
            title(['Surface Plot of the Temperature at time=',num2str(t),'(sec)'])
            colormap("jet");
            colorbar;
            % exportgraphics(gcf,"cooling_fin3d_unsteady_implicit.gif",'Append',true);
            movegui(f2,'west')
        end

        close(progress)
        disp('Solved unsteady state problem using Implicit Theta Method')
end




