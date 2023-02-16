function [X, Y] = setUpMesh(M, l, formfunction)
[dimY, dimX] = size(M);
x = repmat(linspace(0,l,dimX),dimY,1);
y_end = formfunction(x(1,:)/l);
y = zeros(dimY,dimX);
for i =dimX:-1:1
    y(:,i) = linspace(0,y_end(i),dimY);
end

X = x;
Y = flip(y,1);
end