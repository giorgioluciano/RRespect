function [taui, gi] =taugcal(p, Gp1, r, taumin, N)
taui = zeros(1, N);
H = zeros(1, N);
gi = zeros(1,N);
taui(1) = taumin;
for i=1:N-1
    taui(i+1) = taui(i)*r;
end
G1 = 2*Gp1*gamma(1+p)/p/pi*sin(p*pi/2);
for i = 1:N
    H(i) = G1*taui(i)^(-p)/gamma(p);
    gi(i) = H(i)*(r^(p/2) - r^(-p/2))/p;
end
end
    


