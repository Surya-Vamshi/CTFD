function [x, status] = solver_conjugate_gradient_squared(A,b,x_init,max_iter,tol)
status = 0;
x0 = x_init;
x = x0;
r0 = b - A*x0;
u_k = 0;
w_k = 0;
r0_hat = r0;
alpha = 1;
sigma = 1;
for iter = 1:max_iter
    rhok = dot(r0_hat,r0);
    beta = (-1/alpha)*(rhok/sigma);
    v_k = r0 - beta*u_k;
    w_k = v_k - beta*(u_k - beta*w_k);
    c = A*w_k;
    sigma = dot(c, r0_hat);
    alpha = rhok/sigma;
    u_k = v_k - alpha*c;
    x = x + alpha*(v_k + u_k);
    if norm(b-A*x) < tol % Check for convergence
        fprintf('Conjugate Gradient Squared solver converged after %i iterations\n',iter)
        status = 1;
        break;
    end
    r0 = r0 - alpha*A*(v_k + u_k);
end
end