function [y2] = interpolate(x1,y1,x2,x3,y3)
y2 = (((x2-x1)*(y3-y1))/(x3-x1)) + y1;

end