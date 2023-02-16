function x = multigrid_prolongate(xc, dimX, dimY)
% Interpolate xc to the fine grid using a simple injection operator
xc = reshape(xc, dimX, dimY)';
x = zeros(dimX*2, dimY*2);

for i = 1:dimX - 1
    for j = 1:dimY - 1
        x(2 * i, 2 * j) = xc(i, j);
        x(2 * i + 1, 2 * j) = 0.5 * (xc(i, j) + xc(i + 1, j));
        x(2 * i, 2 * j + 1) = 0.5 * (xc(i, j) + xc(i, j + 1));
        x(2 * i + 1, 2 * j + 1) = 0.25 * (xc(i, j) + xc(i + 1, j) + xc(i, j + 1) + xc(i + 1, j + 1));
    end
end

for i = 2:dimX - 1
    x(2 * i - 1, 1) = xc(i, 1);
    x(2 * i, 1) = xc(i, 1);
    x(2 * i - 1, dimY*2) = xc(i, dimY);
    x(2 * i, dimY*2) = xc(i, dimY);
end

for j = 2:dimY - 1
    x(1, 2 * j - 1) = xc(1, j);
    x(1, 2 * j) = xc(1, j);
    x(dimX*2, 2 * j - 1) = xc(dimX, j);
    x(dimX*2, 2 * j) = xc(dimX, j);
end

x = x(:);
end