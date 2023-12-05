function [rho] = marsDrag(altitude)

if altitude > 1000
    rho = 0;
else

    xData = [0.00651138001957105;0.00171986958638280;0.000544930489881507;0.000120080309514958;...
        5.14862489268230e-05;1.20502848981919e-05;1.20996088624937e-06;1.74586938158406e-07;...
        3.62284622377271e-08;6.26707561273697e-09;7.53989100327909e-10;1.08752642421701e-10;...
        1.56681444785521e-11;2.39893891627418e-12;5.95894381854078e-13;2.12668397813016e-13;...
        5.27360453613484e-14;1.99749299976395e-14;6.69652601953330e-15;3.42915378455658e-15;...
        1.98209108523692e-15;1.14523493491467e-15;8.94933835984164e-16;6.58559062363830e-16;...
        5.14723390998580e-16;4.02226097662187e-16;2.95593046396505e-16;2.17353521861868e-16];

    yData = [6.62460567823337;21.7665615141957;29.3375394321765;41.9558359621451;47.0031545741324;...
        62.1451104100945;74.7634069400631;89.9053627760251;97.4763406940062;112.618296529968;...
        132.807570977918;152.996845425867;188.328075709779;218.611987381703;246.372239747634;...
        279.179810725552;329.652996845426;375.078864353312;435.646687697161;491.167192429022;...
        544.164037854889;602.208201892744;665.299684542587;723.343848580442;783.911671924290;...
        847.003154574132;922.712933753943;990.851735015773];

    for i = 1:length(yData)
        if yData(i) > altitude
            break
        end
    end
    x1 = xData(i-1);
    y1 = yData(i-1);
    y2 = altitude;
    x3 = xData(i);
    y3 = yData(i);

    rho = interpolate(y1,x1,y2,y3,x3);

end
end