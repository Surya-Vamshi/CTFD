function [x, status] = solver_steepest_descent(A,b,x0,max_iter,tol)
% Initial guess for iterative solvers
x = x0;

status = 0;
% Main solver loop
for iter = 1:max_iter
    % Update x using the steepest descent method
    r = b - A*x; % Compute the residual
    alpha = (r'*r) / (r'*A*r); % Compute the step size
    x = x + alpha*r; % Update x
    residual = norm(b - A*x); % Compute the new residual
    if residual < tol
        fprintf('Steepest descent solver converged after %i iterations\n',iter)
        status = 1;
        break;
    end
end
if iter == max_iter
    status = 0;
    fprintf('Steepest descent solver did not converged, residual = %f \n',residual)
end
end