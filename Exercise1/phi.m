function [y] = phi(x, L)
%PHI Actual phi function

y = piecewise((x>0)&(x<=1),x,(x>1)&(x<=2),2-x,(x>2)&(x<L),0);
end
