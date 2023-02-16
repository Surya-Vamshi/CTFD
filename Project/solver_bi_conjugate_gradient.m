function [x, status] = solver_bi_conjugate_gradient(A,b,x_init,max_iter,tol)
status = 0;
x0 = x_init;
r0 = b - A*x0;
r0_hat = r0;
rho0 = 1;
alpha = 1;
omega0 = 1;
v0 = zeros(size(x0));
p0 = zeros(size(x0));
for iter = 1:max_iter
    rho1 = dot(r0_hat,r0);
    beta = (rho1/rho0)*(alpha/omega0);
    p1 = r0 + beta*(p0 - omega0*v0);
    v1 = A*p1;
    alpha = rho1/dot(r0_hat, v1);
    h = x0 + alpha*p1;
    if norm(b-A*h) < tol % Check for convergence
        fprintf('Bi Conjugate Gradient solver converged after %i iterations\n',iter - 0.5)
        x = h;
        status = 1;
        break;
    end
    s = r0 - alpha*v1;
    t = A * s;
    omega1 = dot(t,s)/dot(t,t);
    x1 = h + omega1*s;
    if norm(b-A*x1) < tol % Check for convergence
        fprintf('Bi Conjugate Gradient solver converged after %i iterations\n',iter)
        x = x1;
        status = 1;
        break;
    end
    r1 = s - omega1*t;
    rho0 = rho1;
    p0 = p1;
    omega0 = omega1;
    v0 = v1;
    x0 = x1;
    r0 = r1;
end
end

% function [x, status] = solver_bi_conjugate_gradient(A,b,x0,max_iter,tol)
% status = 0;
% % Initial guess for iterative solvers
% x = x0;
% r = b - A*x; % Compute the initial residual
% rho = 1;
% alpha = 1;
% omega = 1;
% p = zeros(size(x));
% v = zeros(size(x));
% % Main solver loop
% for iter = 1:max_iter
%     rho_prev = rho;
%     rho = r'*r;
%     beta = (rho / rho_prev) * (alpha / omega);
%     p = r + beta * (p - omega * v);
%     v = A*p;
%     alpha = rho / (r'*v);
%     s = r - alpha * v;
%     if norm(s) < tol % Check for convergence
%         break;
%     end
%     t = A*s;
%     omega = (t'*s) / (t'*t);
%     x = x + alpha * p + omega * s;
%     r = s - omega * t;
%     if norm(r) < tol
%         fprintf('Bi Conjugate Gradient solver converged after %i iterations\n',iter)
%         status = 1;
%         break;
%     end
% end
% 
% if iter == max_iter
%     status = 0;
%     fprintf('Bi Conjugate Gradient solver did not converged, residual = %f \n',norm(r))
% end
% end

% function [x, status] = solver_bi_conjugate_gradient(A,b,x0,max_iter,tol)
% status = 0;
% Initial guess for iterative solvers
% x = x0;
% r = b - A*x; % Compute the initial residual
% rho = 1;
% alpha = 1;
% omega = 1;
% p = zeros(size(x));
% v = zeros(size(x));
% Main solver loop
% for iter = 1:max_iter
%     for i = 1:length(b)
%         rho_prev = rho;
%         rho = r'*r;
%         beta = (rho / rho_prev) * (alpha / omega);
%         p = r + beta * (p - omega * v);
%         v = A*p;
%         alpha = rho / (r'*v);
%         s = r - alpha * v;
%         if norm(s) < tol % Check for convergence
%             break;
%         end
%         t = A*s;
%         omega = (t'*s) / (t'*t);
%         x = x + alpha * p + omega * s;
%         r = s - omega * t;
%     end
%     if norm(r) < tol
%         fprintf('Bi Conjugate Gradient solver converged after %i iterations\n',iter)
%         status = 1;
%         break;
%     end
% end
% 
% if iter == max_iter
%     status = 0;
%     fprintf('Bi Conjugate Gradient solver did not converged, residual = %f \n',norm(r))
% end
% end
