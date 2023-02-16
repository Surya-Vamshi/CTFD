function [x, status] = solver_gmres( A, b, x_init, max_iter, tol, restart)
x = x_init;
r = b - A * x;
residuals = [];
residuals(1) = norm(r);
status = 0;

for outer_iter = 1:max_iter
    V = zeros(length(b), restart+1);
    H = zeros(restart+1, restart);
    V(:,1) = r / residuals(end);
    for inner_iter = 1:restart
        w = A * V(:,inner_iter);
        for j = 1:inner_iter
            H(j, inner_iter) = w' * V(:,j);
            w = w - H(j, inner_iter) * V(:,j);
        end
        H(inner_iter+1, inner_iter) = norm(w);
        if H(inner_iter+1, inner_iter) < 1e-15
            break;
        end
        V(:,inner_iter+1) = w / H(inner_iter+1, inner_iter);
    end
    y = H(1:inner_iter, 1:inner_iter) \ (residuals(end) * eye(inner_iter,1));
    x = x + V(:,1:inner_iter) * y;
    r = b - A * x;
    residuals(end+1) = norm(r);
    if residuals(end) < tol
        fprintf('Generalized minimal residual solver with restart converged after %i iterations\n',outer_iter);
        status = 1;
        break;
    end 
end
