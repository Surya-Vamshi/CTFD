function [x, status] = solver_gauss_seidel(A,b,x0,max_iter,tol)
% Initial guess for iterative solvers
x = x0;

status = 0;
% Main solver loop
for iter = 1:max_iter
    for j = 1:length(b)
        x(j) = (1/A(j,j))*(b(j) - A(j,:)*x + A(j,j)*x(j));
    end
    eps= norm(b-(A*x))/norm(b);
    if eps < tol
        status = 1;
        fprintf('Gauss Seidel solver converged after %i iterations\n',iter)
        break;
    end
end
if iter == max_iter
    status = 0;
    fprintf('Gauss Seidel  solver did not converged, residual = %f \n',eps)
end
end
