function x = multigrid_wcycle_recursive(A, b, x, R, P, max_level, curr_level)
omega = 2/3;
if curr_level < max_level
    %  Pre-smoothing
    %  x = smoother_gauss_seidel(Ah,b,x,25);
    x = smoother_weightedJacobi(A{curr_level}, b, x, omega, 25);

    % Compute residual
    rh = b - A{curr_level}*x;

    % Restrict residual to coarse grid
    r2h = R{curr_level} * rh;
    e2h = zeros(size(r2h));

    % Recursive call to the W-cycle algorithm
    e2h = multigrid_wcycle_recursive(A, r2h, e2h, R, P, max_level, curr_level+1);

    % Prolongate and correct solution
    x = x + P{curr_level} * e2h;

    %  Re-smoothing
    %  x = smoother_gauss_seidel(Ah,b,x,25);
    x = smoother_weightedJacobi(A{curr_level}, b, x, omega, 25);

    % Compute residual
    rh = b - A{curr_level}*x;

    % Restrict residual to coarse grid
    r2h = R{curr_level} * rh;

    % Recursive call to the W-cycle algorithm
    e2h = multigrid_wcycle_recursive(A, r2h, e2h, R, P, max_level, curr_level+1);

    % Prolongate and correct solution
    x = x + P{curr_level} * e2h;

    %  Post-smoothing
    %  x = smoother_gauss_seidel(Ah,b,x,25);
    x = smoother_weightedJacobi(A{curr_level}, b, x, omega, 25);
else
    %     [x, status] = solver_gauss_seidel(Ah,b,x0,25000,tol);
    x = A{curr_level} \ b;
end
end