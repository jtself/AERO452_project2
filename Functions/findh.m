function h = findh(a, mu, ecc)
% given magnitude r, mu, eccentricity, and true anomaly, this function
% finds the magnitude of h, the specific angular momentum of ANY orbit

h = sqrt(a*mu*(1 - ecc^2));
end