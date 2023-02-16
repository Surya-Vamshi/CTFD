function x = smoother_gauss_seidel(A,b,x0,numIterations)
% Initial guess for iterative solvers
x = x0;

% Main solver loop
for iter = 1:numIterations
    for j = 1:length(b)
        x(j) = (1/A(j,j))*(b(j) - A(j,:)*x + A(j,j)*x(j));
    end
end
end

