clear all
% CODI FINAL!!!

%% DATOS
M =  200;  % número de paneles
f =  0.02;  % curvatura máxima 
p =  0.4;  % posición curvatura máxima
c = 1; %longitud cuerda
xh =  1.0;  % posición hinge 
eta =  0.0; % ángulo deflexión flap
U_inf = 1; % velocidad de corriente libre
x_ref = c/4; % referencia: c/4

%% 2.1 TABLA
[Tabla, CL_alfa_taula, alfa_l0_taula] = Validacion_sin_flap(M, f, p, c, xh, eta, U_inf, x_ref);
