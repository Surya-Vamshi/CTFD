function P = multigrid_prolongate_matrix(dimX_2, dimY_2)
%% Constructing the restriction matrix
dimX = dimX_2/2;
dimY = dimY_2/2;
index = @(ii, jj) ii + (jj-1)*dimX;
index_p = @(ii, jj) ii + (jj-1)*dimX*2;

P = zeros(dimX*dimY*4, dimX*dimY);

for i = 2:dimX - 1
    for j = 2:dimY - 1
        P(index_p(2*i,2*j), index(i,j)) = 1;
        P(index_p(2*i+1, 2*j), index(i,j)) = 0.5;
        P(index_p(2*i-1, 2*j), index(i,j)) = 0.5;
        P(index_p(2*i, 2*j+1), index(i,j)) = 0.5;
        P(index_p(2*i, 2*j-1), index(i,j)) = 0.5;
        P(index_p(2*i-1, 2*j-1), index(i,j)) = 0.25;
        P(index_p(2*i-1, 2*j+1), index(i,j)) = 0.25;
        P(index_p(2*i+1, 2*j+1), index(i,j)) = 0.25;
        P(index_p(2*i+1, 2*j-1), index(i,j)) = 0.25;
    end
end

for i = 1:dimX 
    P(index_p(2*i,1), index(i,1)) = 1;
    P(index_p(2*i-1,1), index(i,1)) = 1;

    P(index_p(2*i,dimY*2), index(i,dimY)) = 1;
    P(index_p(2*i-1,dimY*2), index(i,dimY)) = 1;
end

for j = 1:dimY
    P(index_p(1,2*j), index(1,j)) = 1;
    P(index(1,2*j-1), index(1,j)) = 1;

    P(index_p(dimX*2,2*j), index(dimX,j)) = 1;
    P(index_p(dimX*2,2*j-1), index(dimX,j)) = 1;
end
end