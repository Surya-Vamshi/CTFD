function R = multigrid_restrict_matrix(dimX, dimY)
%% Constructing the restriction matrix
index_r = @(ii, jj) ii + (jj-1)*dimX/2;
index = @(ii, jj) ii + (jj-1)*dimX;

R = zeros(dimX*dimY/4, dimX*dimY);

for i = 2:dimX/2 - 1
    for j = 2:dimY/2 - 1
        R(index_r(i,j), index(2*i,2*j)) = 0.25;
        R(index_r(i,j), index(2*i+1, 2*j)) = 0.125;
        R(index_r(i,j), index(2*i-1, 2*j)) = 0.125;
        R(index_r(i,j), index(2*i, 2*j+1)) = 0.125;
        R(index_r(i,j), index(2*i, 2*j-1)) = 0.125;
        R(index_r(i,j), index(2*i-1, 2*j-1)) = 0.0625;
        R(index_r(i,j), index(2*i-1, 2*j+1)) = 0.0625;
        R(index_r(i,j), index(2*i+1, 2*j+1)) = 0.0625;
        R(index_r(i,j), index(2*i+1, 2*j-1)) = 0.0625;
    end
end

for i = 2:dimX/2 - 1
    R(index_r(i,1), index(2*i,1)) = 0.5;
    R(index_r(i,1), index(2*i-1,1)) = 0.25;
    R(index_r(i,1), index(2*i+1,1)) = 0.25;

    R(index_r(i,dimY/2), index(2*i,dimY)) = 0.5;
    R(index_r(i,dimY/2), index(2*i-1,dimY)) = 0.25;
    R(index_r(i,dimY/2), index(2*i+1,dimY)) = 0.25;  
end

for j = 2:dimY/2 - 1
    R(index_r(1,j), index(1,2*j)) = 0.5;
    R(index_r(1,j), index(1,2*j-1)) = 0.25;
    R(index_r(1,j), index(1,2*j+1)) = 0.25;

    R(index_r(dimX/2,j), index(dimX,2*j)) = 0.5;
    R(index_r(dimX/2,j), index(dimX,2*j-1)) = 0.25;
    R(index_r(dimX/2,j), index(dimX,2*j+1)) = 0.25; 
end


R(index_r(1,1), index(1,1)) = 0.5;
R(index_r(1,1), index(1,2)) = 0.25;
R(index_r(1,1), index(2,1)) = 0.25;

R(index_r(dimX/2,1), index(dimX,1)) = 0.5;
R(index_r(dimX/2,1), index(dimX-1,1)) = 0.25;
R(index_r(dimX/2,1), index(dimX,2)) = 0.25;

R(index_r(1,dimY/2), index(1,dimY)) = 0.5;
R(index_r(1,dimY/2), index(2,dimY)) = 0.25;
R(index_r(1,dimY/2), index(1,dimY-1)) = 0.25;

R(index_r(dimX/2,dimY/2), index(dimX,dimY)) = 0.5;
R(index_r(dimX/2,dimY/2), index(dimX-1,dimY)) = 0.25;
R(index_r(dimX/2,dimY/2), index(dimX,dimY-1)) = 0.25;

end