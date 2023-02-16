function [x, status] = multigrid_2level(Ah, b, x0, dimX, dimY, boundary, max_iter, tol)
% Initialization
x = x0;
status = 0;
omega = 2/3;
iter=0;
restrictor = multigrid_restrict_matrix(dimX, dimY);
prolongator = multigrid_prolongate_matrix(dimX, dimY);

% Main Loop
while(status == 0 && iter<max_iter)
    %  Pre-smoothing
    %  x = smoother_gauss_seidel(Ah,b,x,25);
    x = smoother_weightedJacobi(Ah, b, x, omega, 25);

    % Calculating Residuals
    rh = b - Ah*x;

    % Restrict residual to coarse grid
    % r2h = multigrid_restrict(rh, dimX, dimY);
    r2h = restrictor*rh;

    % Getting Restricted A matrix
    [A2h, ~] = construct_matrix(dimX/2, dimY/2, boundary);
    
    % Solving the course equation using exact method
    e2h = A2h \ r2h;

    % Prolongate and correct solution
    % x = x + multigrid_prolongate(e2h, dimX/2, dimY/2);
    x = x + prolongator*e2h;

    %  Post-smoothing
    x = smoother_weightedJacobi(Ah,b,x,omega,25);

    % Check convergence
    r = b - Ah*x;
    if norm(r) < tol
        fprintf('Multigird V cycle with 2 levels converged after %i iterations\n',iter)
        status = 1;
    end
    iter = iter + 1;
end
if iter == max_iter
    status = 0;
    fprintf('Multigird V cycle with 2 levels did not converge even after %i iterations\n',iter)
    fprintf('Residual = %f \n',norm(r))
end
end
