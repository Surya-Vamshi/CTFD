function [h, q] = arnoldi(A, Q, k)
% Arnoldi Function for Generalized minimal residual method
q = A*Q(:,k);
for i = 1:k
    h(i) = q' * Q(:, i);
    q = q - h(i) * Q(:, i);
end
h(k + 1) = norm(q);
q = q / h(k + 1);
end