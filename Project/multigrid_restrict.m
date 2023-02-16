function xc = multigrid_restrict(x, dimX, dimY)
x = reshape(x, dimX, dimY)';
xc = zeros(dimX/2, dimY/2);
for i = 2:dimX/2 - 1
    for j = 2:dimY/2 - 1
        xc(i,j) = 0.25 * x(2*i, 2*j) + 0.125*(x(2*i-1,2*j)+x(2*i+1,2*j) + x(2*i,2*j-1) + x(2*i,2*j+1)) + 0.0625*(x(2*i-1, 2*j-1) + x(2*i-1,2*j+1) + x(2*i+1,2*j-1) + x(2*i+1,2*j+1));
    end
end

for i = 2:dimX/2 - 1
    xc(i,1) = 0.5*x(2*i, 1) + 0.25*(x(2*i-1, 1) + x(2*i+1, 1));
    xc(i,dimY/2) = 0.5*x(2*i, dimY) + 0.25*(x(2*i-1, dimY) + x(2*i+1, dimY));
end

for j = 2:dimY/2 - 1
    xc(1,j) = 0.5*x(1,2*j) + 0.25*(x(1, 2*j-1) + x(1, 2*j+1));
    xc(dimX/2,j) = 0.5*x(dimX,2*j) + 0.25*(x(dimX, 2*j-1) + x(dimX, 2*j+1));
end

xc = xc(:);
end