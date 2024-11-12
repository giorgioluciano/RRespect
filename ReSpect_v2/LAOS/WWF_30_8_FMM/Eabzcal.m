function [Eabz] =Eabzcal(a, b, z, k)

Eabz = z.^k ./gamma(a .*k + b);

end
    


