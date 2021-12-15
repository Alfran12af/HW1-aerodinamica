function [CLTAT, CMLETAT] = CoeficientesTAT(p, f, c, alfa)


dxdz0 = @(theta) (2*p-1+cos(theta));
dxdz1 = @(theta) (2*p-1+cos(theta)).*cos(theta);
dxdz2 = @(theta) (2*p-1+cos(theta)).*cos(2*theta);


theta_p = acos(1-2*p);
A0 = alfa - 1/pi*(c*f/(p^2)*integral(dxdz0,0,theta_p) + c*f/((1-p)^2)*integral(dxdz0,theta_p, pi));
A1 = 2/pi*(c*f/p^2*integral(dxdz1,0,theta_p) + c*f/((1-p)^2)*integral(dxdz1,theta_p, pi));
A2 = 2/pi*(c*f/p^2*integral(dxdz2,0,theta_p) + c*f/((1-p)^2)*integral(dxdz2,theta_p, pi));

CLTAT = (2*A0+A1)*pi;
CMLETAT = -pi/2*(A0+A1-A2/2);