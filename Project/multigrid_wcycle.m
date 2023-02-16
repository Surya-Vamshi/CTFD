function [x, status] = multigrid_wcycle(Ah, b, x0, dimX, dimY, boundary, max_iter, tol, max_level)
% Initialization
x = x0;
status = 0;
iter=0;
R = cell(max_level);
P = cell(max_level);
A = {Ah};
for n = 2:max_level
    [A{n}, ~] = construct_matrix(dimX/(2^(n-1)), dimY/(2^(n-1)), boundary);
end
for n = 1:max_level
    R{n} = multigrid_restrict_matrix(dimX/(2^(n-1)), dimY/2^(n-1));
    P{n} = multigrid_prolongate_matrix(dimX/(2^(n-1)), dimY/2^(n-1));
end

% Main Loop
while(status == 0 && iter<max_iter)
    x = multigrid_wcycle_recursive(A, b, x, R, P, max_level, 1);
    r = b - Ah*x;
    if norm(r) < tol
        fprintf('Multigird W cycle with %i levels converged after %i iterations\n', max_level, iter)
        status = 1;
    end
    iter = iter + 1;
end
if iter == max_iter
    status = 0;
    fprintf('Multigird W cycle with %i levels did not converge even after %i iterations\n', max_level, iter)
    fprintf('Residual = %f \n',norm(r))
end
end