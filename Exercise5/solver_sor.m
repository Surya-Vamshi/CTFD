function [x, status] = solver_sor(A,b,x0,max_iter,tol,omega)
D = diag(diag(A));
L = tril(A,-1);
N = 1/omega*(D+omega*L);
x = x0;
r = b - A*x;

status = 0;
% Main solver loop
for iter = 1:max_iter
    x = x + N\r;
    r = b - A*x;
    eps=norm(r)/norm(b);
    if eps < tol
        status = 1;
        fprintf('SOR solver converged after %i iterations\n',iter)
        break;
    end
end
if iter == max_iter
    status = 0;
    fprintf('SOR solver did not converged, residual = %f \n',eps)
end
end
