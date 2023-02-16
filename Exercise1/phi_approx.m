function [phi_approx] = phi_approx(a_m_array, L, N, v,t)
%PHI_APPROX Approximate calculation of the Phi using fourier description

phi_approx = zeros(2*L + 1, 1);

for i = -L:1:L
    for m = 1:N
        phi_approx(i + L + 1) = phi_approx(i + L + 1) + a_m_array(m)*sin(m*pi*(i - v*t)/L);
    end
end

end

