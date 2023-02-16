function [x, status] = solver_jacobi(A,b,x0,max_iter,tol)
% Initial guess for iterative solvers
x = x0;
% Calculating jacobi iteration matrix
d = diag(A);
J = A./d;
J = J - diag(diag(J));

% Checking criteria for convergance
s = sum(abs(A), 2) - abs(d);
if any(d <= s)
    warning('Matrix A is not diagonal dominant. Hence, solution may not converge.');
end

status = 0;
% Main solver loop
for iter = 1:max_iter
    x_new = J*x + b./d;
    eps= norm(b-(A*x_new))/norm(b);
    if eps < tol
        fprintf('Jacobi solver converged after %i iterations\n',iter)
        status = 1;
        break;
    end
    x = x_new;
end
if iter == max_iter
    status = 0;
    fprintf('Jacobi solver did not converged, residual = %f \n',eps)
end
end

