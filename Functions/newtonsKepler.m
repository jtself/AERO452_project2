function [eccentricAnomaly] = newtonsKepler(Me,ecc)
% Finds the Eccentric Anomaly through Newtons Method
% Must pass in Mean Anomaly and Eccentricity

% Tolerance
TOL = 1e-8;

% Finds initial guess
if Me <= pi
    Eguess = Me + ecc/2;
else
    Eguess = Me - ecc/2;
end

% Builds function F and Fprime
f =  @(E) Me - E + ecc*sin(E);
fprime = @(E) -1 + ecc*cos(E);

% Newtons Method
x_0 = Eguess;
x_1 = x_0 -(f(x_0)/fprime(x_0));
x = [x_0; x_1];
err = abs(x_1 - x_0);

while err > TOL
    x_0 = x_1;
    x_1 = x_0 - f(x_0)/fprime(x_0);
    err = abs(x_1 - x_0);
    x = [x; x_1];
end

eccentricAnomaly = x(end);
end