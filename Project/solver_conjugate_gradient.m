function [x, status] = solver_conjugate_gradient(A,b,x0,max_iter,tol)
% Initial guess for iterative solvers
% [x, ~, relrses, iter] = bicgstab(A,b,tol,max_iter);
% x = A\b;
% [x,flag, relres, iter] = bicgstab(A,b,tol,max_iter,[],[],x0);
% status = ~flag;
% if(status==0)
%     fprintf('Conjugate Gradient solver did not converge, residual = %f \n',relres)
% else
%     fprintf('Conjugate Gradient solver converged after %i iterations\n',iter)
% end
x = x0;
d = A*x - b;
p = d;
iter = 1;
dd = d'*d; % create a variable to store the scalar product

% Main solver loop
while (iter<max_iter && norm(d)>tol)

    % save the only matrix-vector product of the iteration
    ap = A*p;

    % amplitude of the step (1st scalar product for the current iteration)
    alpha = dd/(p'*ap);

    % new approximated solution x
    x(:,iter+1) = x(:,iter) - alpha*p;

    % updating the gradient
    d = d - alpha*ap;

    % definition of beta using the old d and the new one
    % (2nd and last scalar product for the current iteration)
    dd_new = d'*d;
    beta = dd_new/dd;
    dd = dd_new;

    % updating of direction d
    p = d + beta*p;

    iter = iter + 1;
end

iter = iter - 1;
x = x(:,iter);

if iter == max_iter
    status = 0;
    fprintf('Conjugate Gradient solver did not converged, residual = %f \n',norm(r))
else
    fprintf('Conjugate Gradient solver converged after %i iterations\n',iter)
    status = 1;
end
end