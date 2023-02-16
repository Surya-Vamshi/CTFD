function [stencil, b] = stamp(i, j, X, Y, b, alpha, TD, Tinf, lambda, q_dot_sym, boundary, verbose)
%stecil calculate the linear equation for node (i, j)

%  input:
%    i         node number in x direction
%    j         node number in y direction
%    X         x position of the nodes
%    Y         y position of the nodes
%    b         right-hand side value for node (i,j)
%    alpha     alpha
%    Tinf      Tinf for Robin BC
%    boundary  defines the boundary conditions
%    verbose   verbositiy level
%
%  output:
%    stecil     linear equation for node (i,j)
%    b         new right-hand side value for node (i,j)


% Init

n = size(X, 1);
m = size(X, 2);
stencil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

% Determine the node positon
if i==1 && j==1
    nodePosition = 'NW';
elseif i==n && j==m
    nodePosition = 'SE';
elseif i==1 && j==m
    nodePosition = 'NE';
elseif i==n && j==1
    nodePosition = 'SW';
elseif j == 1
    nodePosition = 'West';
elseif j == size(X,2)
    nodePosition = 'East';
elseif i == size(X,1)
    nodePosition = 'South';
elseif i == 1
    nodePosition = 'North';
else
    nodePosition = 'inner Node';
end


% Calculate the equation for the correct node position
switch nodePosition

    case 'inner Node'


        % Nomenclature:
        %
        %    NW(i-1,j-1)   Nw -  N(i-1,j) -  Ne     NE(i-1,j+1)
        %
        %                 |                 |
        %
        %       nW - - - - nw ------ n ------ ne - - - nE
        %                 |                 |
        %       |         |        |        |       |
        %                 |                 |
        %   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)
        %                 |                 |
        %       |         |        |        |       |
        %                 |                 |
        %      sW - - - - sw ------ s ------ se - - - sE
        %
        %                 |                 |
        %
        %   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)
        %
        % Indexing of stecil:

        %    D_4 - D_1 - D2
        %     |     |     |
        %    D_3 - D_0 - D3
        %     |     |     |
        %    D_2 -  D1 - D4

        % Principal node coordinates
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1);
        y_N  = Y(i-1,j);     x_N  = X(i-1,j);
        y_W  = Y(i,j-1);     x_W  = X(i,j-1);
        y_E  = Y(i,j+1);     x_E  = X(i,j+1);
        y_S  = Y(i+1,j);     x_S  = X(i+1,j);
        y_P  = Y(i,j);       x_P  = X(i,j);
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1);
        y_SE = Y(i+1,j+1);   x_SE = X(i+1,j+1);
        y_NE = Y(i-1,j+1);   x_NE = X(i-1,j+1);

        % Auxiliary node coordinates
        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;
        y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
        y_nW = (y_NW + y_W)/2;  x_nW = (x_NW + x_W)/2;
        y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
        y_Sw = (y_SW + y_S)/2;  x_Sw = (x_SW + x_S)/2;
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
        y_sW = (y_SW + y_W)/2;  x_sW = (x_SW + x_W)/2;
        y_sE = (y_SE + y_E)/2;  x_sE = (x_SE + x_E)/2;

        y_n  = (y_N + y_P)/2;   x_n  = (x_N + x_P)/2;
        y_s  = (y_S + y_P)/2;   x_s  = (x_S + x_P)/2;
        y_e  = (y_E + y_P)/2;   x_e  = (x_E + x_P)/2;
        y_w  = (y_W + y_P)/2;   x_w  = (x_W + x_P)/2;

        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
        y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
        y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;

        % Inter node distances
        % Around s
        dy_Sw_Se = y_Se - y_Sw;   dx_Sw_Se = x_Se - x_Sw;
        dy_Se_e  = y_e  - y_Se;   dx_Se_e  = x_e  - x_Se;
        dy_e_w   = y_w  - y_e;    dx_e_w   = x_w  - x_e;
        dy_w_Sw  = y_Sw - y_w;    dx_w_Sw  = x_Sw - x_w;

        % Around e
        dy_s_sE  = y_sE - y_s;    dx_s_sE  = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;   dx_sE_nE = x_nE - x_sE;
        dy_nE_n  = y_n  - y_nE;   dx_nE_n  = x_n  - x_nE;
        dy_n_s   = y_s  - y_n;    dx_n_s   = x_s  - x_n;

        % Around n
        dy_w_e   = y_e  - y_w;    dx_w_e   = x_e  - x_w;
        dy_e_Ne  = y_Ne - y_e;    dx_e_Ne  = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;   dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w  = y_w  - y_Nw;   dx_Nw_w  = x_w  - x_Nw;

        % Around w
        dy_sW_s  = y_s  - y_sW;   dx_sW_s  = x_s  - x_sW;
        dy_s_n   = y_n  - y_s;    dx_s_n   = x_n  - x_s;
        dy_n_nW  = y_nW - y_n;    dx_n_nW  = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;   dx_nW_sW = x_sW - x_nW;

        % Around P
        dy_sw_se = y_se - y_sw;   dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se;   dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw;   dx_nw_sw = x_sw - x_nw;

        % Areas

        S_P = 0.5*abs(x_ne * y_se - x_se * y_ne + x_se * y_sw - x_sw * y_se...
            + x_sw * y_nw - x_nw * y_sw + x_nw * y_ne - x_ne * y_nw);

        S_n  = 0.5*abs(x_Ne * y_e - x_e * y_Ne + x_e * y_w - x_w * y_e...
            +  x_w * y_Nw - x_Nw * y_w + x_Nw * y_Ne - x_Ne * y_Nw);

        S_e  = 0.5*abs(x_nE * y_sE - x_sE * y_nE + x_sE * y_s - x_s * y_sE...
            +  x_s * y_n - x_n * y_s + x_n * y_nE - x_nE * y_n);

        S_w  = 0.5*abs(x_n * y_s - x_s * y_n + x_s * y_sW - x_sW * y_s...
            +  x_sW * y_nW - x_nW * y_sW + x_nW * y_n - x_n * y_nW);

        S_s = 0.5*abs(x_e * y_Se - x_Se * y_e + x_Se * y_Sw - x_Sw * y_Se...
            + x_Sw * y_w - x_w * y_Sw + x_w * y_e - x_e * y_w);



        %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$

        build_inner

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        % P
        stencil(index(i, j)) = D0;
        % East
        stencil(index(i,j+1)) = D3;
        % West
        stencil(index(i,j-1)) = D_3;
        % South
        stencil(index(i+1,j)) = D1;
        % North
        stencil(index(i-1,j)) = D_1;
        % NW
        stencil(index(i-1,j-1)) = D_4;
        % NE
        stencil(index(i-1,j+1)) = D2;
        % SW
        stencil(index(i+1,j-1)) = D_2;
        % SE
        stencil(index(i+1,j+1)) = D4;


    case 'South'
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1);
        y_N  = Y(i-1,j);     x_N  = X(i-1,j);
        y_NE = Y(i-1,j+1);   x_NE = X(i-1,j+1);
        y_W  = Y(i,j-1);     x_W  = X(i,j-1);
        y_P  = Y(i,j);       x_P  = X(i,j);
        y_E  = Y(i,j+1);     x_E  = X(i,j+1);

        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;
        y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
        y_nW = (y_NW + y_W)/2;  x_nW = (x_NW + x_W)/2;
        y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
        y_n  = (y_N + y_P)/2;   x_n  = (x_N + x_P)/2;
        y_e  = (y_E + y_P)/2;   x_e  = (x_E + x_P)/2;
        y_w  = (y_W + y_P)/2;   x_w  = (x_W + x_P)/2;
        y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;

        dy_w_e   = y_e  - y_w;      dx_w_e   = x_e  - x_w;
        dy_e_ne  = y_ne - y_e;      dx_e_ne  = x_ne - x_e;
        dy_ne_nw = y_nw - y_ne;     dx_ne_nw = x_nw - x_ne;
        dy_nw_w  = y_w  - y_nw;     dx_nw_w  = x_w  - x_nw;
        dy_P_E   = y_E  - y_P;      dx_P_E   = x_E  - x_P;
        dy_E_nE  = y_nE - y_E;      dx_E_nE  = x_nE - x_E;
        dy_nE_n  = y_n  - y_nE;     dx_nE_n  = x_n  - x_nE;
        dy_n_P   = y_P  - y_n;      dx_n_P   = x_P  - x_n;
        dy_e_Ne  = y_Ne - y_e;      dx_e_Ne  = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;     dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w  = y_w  - y_Nw;     dx_Nw_w  = x_w  - x_Nw;
        dy_W_P   = y_P  - y_W;      dx_W_P   = x_P  - x_W;
        dy_P_n   = y_n  - y_P;      dx_P_n   = x_n  - x_P;
        dy_n_nW  = y_nW - y_n;      dx_n_nW  = x_nW - x_n;
        dy_nW_W  = y_W  - y_nW;     dx_nW_W  = x_W  - x_nW;

        S_v = 0.5*abs(x_ne*y_e - y_ne*x_e + x_e*y_w - x_w*y_e ...
            + x_w*y_nw - x_nw*y_w + x_nw*y_ne - x_ne*y_nw);

        S_vw = 0.5*abs(x_n*y_P - x_P*y_n + x_P*y_W - x_W*y_P...
            + x_W*y_nW - x_nW*y_W + x_nW*y_n - x_n*y_nW);

        S_ve = 0.5*abs(x_n*y_nE - x_nE*y_n + x_nE*y_E - x_E*y_nE...
            + x_E*y_P - x_P*y_E + x_P*y_n - x_n*y_P);

        S_n  = 0.5*abs(x_Ne * y_e - x_e * y_Ne + x_e * y_w - x_w * y_e...
            +  x_w * y_Nw - x_Nw * y_w + x_Nw * y_Ne - x_Ne * y_Nw);

        % P
        D0=((dx_e_ne*(dx_P_E/2 + (3*dx_n_P)/4 + dx_nE_n/4))/S_ve + (dx_nw_w*(dx_W_P/2 + (3*dx_P_n)/4 + dx_n_nW/4))/S_vw + (dy_e_ne*(dy_P_E/2 + (3*dy_n_P)/4 + dy_nE_n/4))/S_ve + (dy_nw_w*(dy_W_P/2 + (3*dy_P_n)/4 + dy_n_nW/4))/S_vw + (dx_ne_nw*(dx_w_e + dx_e_Ne/4 + dx_Nw_w/4))/S_n + (dy_ne_nw*(dy_w_e + dy_e_Ne/4 + dy_Nw_w/4))/S_n)/S_v;

        % East
        D3=((dx_e_ne*(dx_P_E/2 + (3*dx_E_nE)/4 + dx_nE_n/4))/S_ve + (dy_e_ne*(dy_P_E/2 + (3*dy_E_nE)/4 + dy_nE_n/4))/S_ve + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_v;

        % NE
        D2=((dx_e_ne*(dx_E_nE/4 + dx_nE_n/4))/S_ve + (dy_e_ne*(dy_E_nE/4 + dy_nE_n/4))/S_ve + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_v;

        % North
        D_1=((dx_e_ne*(dx_n_P/4 + dx_nE_n/4))/S_ve + (dx_nw_w*(dx_P_n/4 + dx_n_nW/4))/S_vw + (dy_e_ne*(dy_n_P/4 + dy_nE_n/4))/S_ve + (dy_nw_w*(dy_P_n/4 + dy_n_nW/4))/S_vw + (dx_ne_nw*(dx_e_Ne/4 + dx_Nw_w/4 + dx_Ne_Nw))/S_n + (dy_ne_nw*(dy_e_Ne/4 + dy_Nw_w/4 + dy_Ne_Nw))/S_n)/S_v;

        % NW
        D_4=((dx_nw_w*(dx_nW_W/4 + dx_n_nW/4))/S_vw + (dy_nw_w*(dy_nW_W/4 + dy_n_nW/4))/S_vw + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dy_Nw_w*dy_ne_nw)/(4*S_n))/S_v;

        % West
        D_3=((dx_nw_w*(dx_W_P/2 + (3*dx_nW_W)/4 + dx_n_nW/4))/S_vw + (dy_nw_w*(dy_W_P/2 + (3*dy_nW_W)/4 + dy_n_nW/4))/S_vw + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dy_Nw_w*dy_ne_nw)/(4*S_n))/S_v;

        stencil(index(i,j)) = D0;
        stencil(index(i,j+1)) = D3;
        stencil(index(i,j-1)) = D_3;
        stencil(index(i-1,j)) = D_1;
        stencil(index(i-1,j-1)) = D_4;
        stencil(index(i-1,j+1)) = D2;
        b = q_dot_sym*sqrt(dx_w_e^2 + dy_w_e^2)/S_v;

        case 'East'
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1);
        y_N  = Y(i-1,j);     x_N  = X(i-1,j);
        y_W  = Y(i,j-1);     x_W  = X(i,j-1);
        y_P  = Y(i,j);       x_P  = X(i,j);
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1);
        y_S  = Y(i+1,j);     x_S  = X(i+1,j);

        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;
        y_nW = (y_NW + y_W)/2;  x_nW = (x_NW + x_W)/2;
        y_sW = (y_SW + y_W)/2;  x_sW = (x_SW + x_W)/2;
        y_Sw = (y_SW + y_S)/2;  x_Sw = (x_SW + x_S)/2;
        y_n  = (y_N + y_P)/2;   x_n  = (x_N + x_P)/2;
        y_w  = (y_W + y_P)/2;   x_w  = (x_W + x_P)/2;
        y_s  = (y_S + y_P)/2;   x_s  = (x_S + x_P)/2;
        y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;
        y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;

        dy_s_n   = y_n  - y_s;      dx_s_n   = x_n  - x_s;
        dy_n_nw  = y_nw - y_n;      dx_n_nw  = x_nw - x_n;
        dy_nw_sw = y_sw - y_nw;     dx_nw_sw = x_sw - x_nw;
        dy_sw_s  = y_s  - y_sw;     dx_sw_s  = x_s  - x_sw;
        dy_P_N   = y_N  - y_P;      dx_P_N   = x_N  - x_P;
        dy_N_Nw  = y_Nw - y_N;      dx_N_Nw  = x_Nw - x_N;
        dy_Nw_w  = y_w  - y_Nw;     dx_Nw_w  = x_w  - x_Nw;
        dy_w_P   = y_P  - y_w;      dx_w_P   = x_P  - x_w;
        dy_n_nW  = y_nW - y_n;      dx_n_nW  = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;     dx_nW_sW = x_sW - x_nW;
        dy_sW_s  = y_s  - y_sW;     dx_sW_s  = x_s  - x_sW;
        dy_S_P   = y_P  - y_S;      dx_S_P   = x_P  - x_S;
        dy_P_w   = y_w  - y_P;      dx_P_w   = x_w  - x_P;
        dy_w_Sw  = y_Sw - y_w;      dx_w_Sw  = x_Sw - x_w;
        dy_Sw_S  = y_S  - y_Sw;     dx_Sw_S  = x_S  - x_Sw;

        S_v = 0.5*abs(x_n*y_s - x_s*y_n + x_s*y_sw - x_sw*y_s...
            + x_sw*y_nw - x_nw*y_sw + x_nw*y_n - x_n*y_nw);

        S_nv = 0.5*abs(x_P*y_w - x_w*y_P + x_w*y_nW - x_nW*y_w...
            + x_nW*y_N - x_N*y_nW + x_N*y_P - x_P*y_N);

        S_sv = 0.5*abs(x_P*y_w - x_w*y_P + x_w*y_sW - x_sW*y_w...
            + x_sW*y_S - x_S*y_sW + x_S*y_P - x_P*y_S);

        S_w  = 0.5*abs(x_n * y_s - x_s * y_n + x_s * y_sW - x_sW * y_s...
            +  x_sW * y_nW - x_nW * y_sW + x_nW * y_n - x_n * y_nW);


        % P
        D0=((dx_n_nw*(dx_P_N/2 + (3*dx_w_P)/4 + dx_Nw_w/4))/S_nv + (dx_sw_s*(dx_S_P/2 + (3*dx_P_w)/4 + dx_w_Sw/4))/S_sv + (dy_n_nw*(dy_P_N/2 + (3*dy_w_P)/4 + dy_Nw_w/4))/S_nv + (dy_sw_s*(dy_S_P/2 + (3*dy_P_w)/4 + dy_w_Sw/4))/S_sv + (dx_nw_sw*(dx_s_n + dx_n_nW/4 + dx_sW_s/4))/S_w + (dy_nw_sw*(dy_s_n + dy_n_nW/4 + dy_sW_s/4))/S_w)/S_v;

        % North
        D_1=((dx_n_nw*(dx_P_N/2 + (3*dx_N_Nw)/4 + dx_Nw_w/4))/S_nv + (dy_n_nw*(dy_P_N/2 + (3*dy_N_Nw)/4 + dy_Nw_w/4))/S_nv + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_v;

        % NW
        D_4=((dx_n_nw*(dx_N_Nw/4 + dx_Nw_w/4))/S_nv + (dy_n_nw*(dy_N_Nw/4 + dy_Nw_w/4))/S_nv + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_v;

        % West
        D_3=((dx_n_nw*(dx_w_P/4 + dx_Nw_w/4))/S_nv + (dx_sw_s*(dx_P_w/4 + dx_w_Sw/4))/S_sv + (dy_n_nw*(dy_w_P/4 + dy_Nw_w/4))/S_nv + (dy_sw_s*(dy_P_w/4 + dy_w_Sw/4))/S_sv + (dx_nw_sw*(dx_n_nW/4 + dx_sW_s/4 + dx_nW_sW))/S_w + (dy_nw_sw*(dy_n_nW/4 + dy_sW_s/4 + dy_nW_sW))/S_w)/S_v;

        % SW
        D_2=((dx_sw_s*(dx_Sw_S/4 + dx_w_Sw/4))/S_sv + (dy_sw_s*(dy_Sw_S/4 + dy_w_Sw/4))/S_sv + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_v;

        % South
        D1=((dx_sw_s*(dx_S_P/2 + (3*dx_Sw_S)/4 + dx_w_Sw/4))/S_sv + (dy_sw_s*(dy_S_P/2 + (3*dy_Sw_S)/4 + dy_w_Sw/4))/S_sv + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_v;

        stencil(index(i,j)) = D0 - alpha*sqrt(dx_s_n^2 + dy_s_n^2)/(lambda*S_v);
        stencil(index(i-1,j)) = D_1;
        stencil(index(i-1,j-1)) = D_4;
        stencil(index(i,j-1)) = D_3;
        stencil(index(i+1,j-1)) = D_2;
        stencil(index(i+1,j)) = D1;
        b = -(Tinf*alpha/lambda)*sqrt(dx_s_n^2 + dy_s_n^2)/S_v;

    case 'North'
        y_W  = Y(i,j-1);     x_W  = X(i,j-1);
        y_P  = Y(i,j);       x_P  = X(i,j);
        y_E  = Y(i,j+1);     x_E  = X(i,j+1);
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1);
        y_S  = Y(i+1,j);     x_S  = X(i+1,j);
        y_SE = Y(i+1,j+1);   x_SE = X(i+1,j+1);

        y_sW = (y_SW + y_W)/2;  x_sW = (x_SW + x_W)/2;
        y_sE = (y_SE + y_E)/2;  x_sE = (x_SE + x_E)/2;
        y_Sw = (y_SW + y_S)/2;  x_Sw = (x_SW + x_S)/2;
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;

        y_e  = (y_E + y_P)/2;   x_e  = (x_E + x_P)/2;
        y_w  = (y_W + y_P)/2;   x_w  = (x_W + x_P)/2;
        y_s  = (y_S + y_P)/2;   x_s  = (x_S + x_P)/2;

        y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;
        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;

        dy_Sw_Se = y_Se - y_Sw;   dx_Sw_Se = x_Se - x_Sw;
        dy_w_Sw  = y_Sw - y_w;    dx_w_Sw  = x_Sw - x_w;
        dy_s_sE  = y_sE - y_s;    dx_s_sE  = x_sE - x_s;
        dy_w_e   = y_e  - y_w;    dx_w_e   = x_e  - x_w;
        dy_sW_s  = y_s  - y_sW;   dx_sW_s  = x_s  - x_sW;
        dy_sw_se = y_se - y_sw;   dx_sw_se = x_se - x_sw;

        dy_e_se  = y_se - y_e;    dx_e_se  = x_se - x_e;
        dy_se_sw = -dy_sw_se;     dx_se_sw = -dx_sw_se;
        dy_sw_w  = y_w - y_sw;    dx_sw_w  = x_w - x_sw;
        dy_P_E   = y_E - y_P;     dx_P_E   = x_E - x_P;
        dy_E_sE  = y_sE - y_E;    dx_E_sE  = x_sE - x_E;
        dy_sE_s  = -dy_s_sE;      dx_sE_s  = -dx_s_sE;
        dy_s_P   = y_P - y_s;     dx_s_P   = x_P - x_s;
        dy_e_Se  = y_Se - y_e;    dx_e_Se  = x_Se - x_e;
        dy_Se_Sw = -dy_Sw_Se;     dx_Se_Sw = -dx_Sw_Se;
        dy_Sw_w = -dy_w_Sw;       dx_Sw_w  = -dx_w_Sw;
        dy_W_P  = y_P - y_W;      dx_W_P   = x_P - x_W;
        dy_P_s  = -dy_s_P;        dx_P_s   = -dx_s_P;
        dy_s_sW = -dy_sW_s;       dx_s_sW  = -dx_sW_s;
        dy_sW_W = y_W - y_sW;     dx_sW_W  = x_W - x_sW;

        S_v = 0.5*abs(x_e*y_se - x_se*y_e + x_se*y_sw - x_sw*y_se...
            + x_sw*y_w - x_w*y_sw + x_w*y_e - x_e*y_w);


        S_vw = 0.5*abs(x_P*y_s - x_s*y_P + x_s*y_sW - x_sW*y_s...
            + x_sW*y_W - x_W*y_sW + x_W*y_P - x_P*y_W);

        S_ve = 0.5*abs(x_P*y_s - x_s*y_P + x_s*y_sE - x_sE*y_s...
            + x_sE*y_E - x_E*y_sE + x_E*y_P - x_P*y_E);

        S_s = 0.5*abs(x_e * y_Se - x_Se * y_e + x_Se * y_Sw - x_Sw * y_Se...
            + x_Sw * y_w - x_w * y_Sw + x_w * y_e - x_e * y_w);

        % P
        D0=((dx_e_se*(dx_P_E/2 + (3*dx_s_P)/4 + dx_sE_s/4))/S_ve + (dx_sw_w*(dx_W_P/2 + (3*dx_P_s)/4 + dx_s_sW/4))/S_vw + (dy_e_se*(dy_P_E/2 + (3*dy_s_P)/4 + dy_sE_s/4))/S_ve + (dy_sw_w*(dy_W_P/2 + (3*dy_P_s)/4 + dy_s_sW/4))/S_vw + (dx_se_sw*(dx_w_e + dx_e_Se/4 + dx_Sw_w/4))/S_s + (dy_se_sw*(dy_w_e + dy_e_Se/4 + dy_Sw_w/4))/S_s)/S_v;

        % East
        D3=((dx_e_se*(dx_P_E/2 + (3*dx_E_sE)/4 + dx_sE_s/4))/S_ve + (dy_e_se*(dy_P_E/2 + (3*dy_E_sE)/4 + dy_sE_s/4))/S_ve + (dx_e_Se*dx_se_sw)/(4*S_s) + (dy_e_Se*dy_se_sw)/(4*S_s))/S_v;

        % SE
        D4=((dx_e_se*(dx_E_sE/4 + dx_sE_s/4))/S_ve + (dy_e_se*(dy_E_sE/4 + dy_sE_s/4))/S_ve + (dx_e_Se*dx_se_sw)/(4*S_s) + (dy_e_Se*dy_se_sw)/(4*S_s))/S_v;

        % South
        D1=((dx_e_se*(dx_s_P/4 + dx_sE_s/4))/S_ve + (dx_sw_w*(dx_P_s/4 + dx_s_sW/4))/S_vw + (dy_e_se*(dy_s_P/4 + dy_sE_s/4))/S_ve + (dy_sw_w*(dy_P_s/4 + dy_s_sW/4))/S_vw + (dx_se_sw*(dx_e_Se/4 + dx_Sw_w/4 + dx_Se_Sw))/S_s + (dy_se_sw*(dy_e_Se/4 + dy_Sw_w/4 + dy_Se_Sw))/S_s)/S_v;

        % SW
        D_2=((dx_sw_w*(dx_sW_W/4 + dx_s_sW/4))/S_vw + (dy_sw_w*(dy_sW_W/4 + dy_s_sW/4))/S_vw + (dx_Sw_w*dx_se_sw)/(4*S_s) + (dy_Sw_w*dy_se_sw)/(4*S_s))/S_v;

        % West
        D_3=((dx_sw_w*(dx_W_P/2 + (3*dx_sW_W)/4 + dx_s_sW/4))/S_vw + (dy_sw_w*(dy_W_P/2 + (3*dy_sW_W)/4 + dy_s_sW/4))/S_vw + (dx_Sw_w*dx_se_sw)/(4*S_s) + (dy_Sw_w*dy_se_sw)/(4*S_s))/S_v;

        stencil(index(i,j)) = D0 - alpha*sqrt(dx_w_e^2 + dy_w_e^2)/(lambda*S_v);
        stencil(index(i+1,j)) = D1;
        stencil(index(i,j-1)) = D_3;
        stencil(index(i,j+1)) = D3;
        stencil(index(i+1,j-1)) = D_2;
        stencil(index(i+1,j+1)) = D4;
        b = -(Tinf*alpha/lambda)*sqrt(dx_w_e^2 + dy_w_e^2)/S_v;

    case 'West'
        stencil(index(i,j)) = 1;

    case 'NE'
        stencil(index(i,j)) = 1;    b = Tinf;

    case 'SE'
        stencil(index(i,j)) = 1;    b = Tinf;

    case 'SW'
        stencil(index(i,j)) = 1;    b = TD.west;

    case 'NW'
        stencil(index(i,j)) = 1;    b = TD.west;

end



