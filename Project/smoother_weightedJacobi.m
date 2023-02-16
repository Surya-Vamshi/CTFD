function x_new = smoother_weightedJacobi(A, b, x, omega, numIterations)
    % Compute the diagonal of A
    D = diag(A);
    % Inverse of the diagonal
    Dinv = 1./D;
    for i = 1:numIterations
        % Compute the residual
        r = b - A*x;
        % Scale the residual by the inverse of the diagonal
        z = Dinv.*r;
        % Update the solution
        x = x + omega*z;
    end
    x_new = x;
end