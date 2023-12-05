function [rho] = marsDrag(altitude)

rho =  0.02*exp(-altitude/11.1);

end