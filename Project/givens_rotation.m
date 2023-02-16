function [h, cs_k, sn_k] = givens_rotation(h, cs, sn, k)
% Applying Givens Rotation for Generalized minimal residual method
for i = 1:k-1
    temp   =  cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
end

t = sqrt(h(k)^2 + h(k + 1)^2);

cs_k = h(k) / t;
sn_k = h(k + 1) / t;

h(k) = cs_k * h(k) + sn_k * h(k + 1);
h(k + 1) = 0.0;
end
