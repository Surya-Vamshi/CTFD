% Author : Vikas Kurapati

clear all;

%% Stage 1
L = 20;
syms x m;

a_m = 2 * int(phi(x, L) * sin(m * pi * x / L), x, [0, L])/L ;

%disp(a_m)

N = 500;

for m = 1:N
    a_m_array(m) = double(subs(a_m, m));
end

%disp(a_m_array)

%% stage 2

phi_approx_ = phi_approx(a_m_array, L, N, 0, 0);

figure(1)
plot(-L:1:L, phi_approx_);
title("Approximation of \phi");
ylabel("\phi");
xlabel("space coordinate");
savefig("Approximation-phi.fig");

%% Stage 3
t_end = 100;
dt = 1;
v = 1;

for t = 0:dt:t_end
    phi_approx_ = phi_approx(a_m_array, L, N, v, t);
    figure(2)
    filename = "convection.gif";
    plot(-L:1:L, phi_approx_);
    title("Convection");
    ylabel("\phi");
    xlabel("space coordinate");
    set(gca, 'YLim', [-1,1]);
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 8);

    if t == 0
        imwrite(imind, cm, filename, 'gif', 'LoopCount', inf);
    else
        imwrite(imind, cm, filename, "gif", 'WriteMode','append', 'DelayTime',0);
    end
end

%% Stage 4
D = 0.1;
t_end = 5;
dt = 0.1;

for t=0:dt:t_end
    phi_diff = zeros(2*L + 1, 1);
    for i = -L:1:L
        for m = 1:N
            phi_diff(i+L+1) = phi_diff(i+L+1) + a_m_array(m)*exp(-D*t*(m*pi/L)^2)*sin(m*pi*i/L);
        end
    end
    
    figure(3)
    filename = "diffusion.gif";
    plot(-L:1:L, phi_diff);
    title("Diffusion");
    ylabel("\phi");
    xlabel("space coordinate");
    set(gca, 'YLim', [-1,1]);
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 8);

    if t==0
        imwrite(imind, cm, filename, 'gif', 'LoopCount',inf);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode','append', 'DelayTime',0);
    end
end

%% Stage 5
D = 0.1;
v = 3;
t_end = 5;
dt = 0.1;

for t = 0:dt:t_end
    phi_convdiff = zeros(2*L + 1, 1);
    for i = -L:1:L
        for m = 1:N
            phi_convdiff(i+L+1) = phi_convdiff(i + L + 1) + a_m_array(m)*exp(-D*t*(m*pi/L)^2)*sin(m*pi*(i-v*t)/L);
        end
    end

    figure(4)
    filename = "convection-diffusion.gif";
    plot(-L:1:L, phi_convdiff);
    title("Convection Diffusion");
    ylabel("\phi");
    xlabel("space coordinate");
    set(gca, 'YLim', [-1,1]);
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 8);

    if t==0
        imwrite(imind, cm, filename, 'gif', 'LoopCount', inf);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode','append', 'DelayTime',0);
    end
end