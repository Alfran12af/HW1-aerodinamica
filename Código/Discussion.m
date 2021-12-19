clear all

%% DATOS
M =  200;  % número de paneles
%f =  0.02;  % curvatura máxima 
%p =  0.4;  % posición curvatura máxima
c = 1; %longitud cuerda

alfa = 4;% ángulo de ataque del análisis
alfa = alfa*pi/180; % conversión de grados a rad

xh =  1.0;  % posición hinge 
eta =  0.0; % ángulo deflexión flap

U_inf = 1; % velocidad de corriente libre
x_ref = 0.25; % referencia: LE

f = 0:0.01:0.06;
p = 0.05:0.05:0.6;
CM0 = zeros(length(f), length(p));
alfa_l0 = zeros(length(f), length(p));

%% C_M0 Y ALFA_L0

for i = 1:length(f) % análisis para cada curvatura máxima
    for j = 1:length(p) % análisis de la posición de cada curvatura máxima
        [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f(i), p(j), c, xh, eta); 
        [G] = Circulacion(M, xcont, xvort, pnorm, alfa, U_inf);
        [CLDVM, CMLEDVM] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa, pchord, c, coord);
        [Tabla, CL_alfa_grados, alfa_l0_asd] = Validacion_sin_flap(M, f(i), p(j), c, xh, eta, U_inf, x_ref);
        
        CM0(i,j) = CMLEDVM; % coeficiente de momento libre calculado por DVM
        alfa_l0(i,j) = alfa_l0_asd; % ángulo de sustentación nula por DVM
    end
end

%% FIGURAS

figure
plot(p(1,:), CM0(1,:), 'LineWidth', 1.3);
hold on
plot(p(1,:), CM0(2,:), 'LineWidth', 1.3);
plot(p(1,:), CM0(3,:), 'LineWidth', 1.3);
plot(p(1,:), CM0(4,:), 'LineWidth', 1.3);
plot(p(1,:), CM0(5,:), 'LineWidth', 1.3);
plot(p(1,:), CM0(6,:), 'LineWidth', 1.3);
plot(p(1,:), CM0(7,:), 'LineWidth', 1.3);
legend({'f=0.00','f=0.01','f=0.02','f=0.03','f=0.04','f=0.05','f=0.06'}, 'Location','southwest');
xlabel('Posición de camber máxima (p)');
ylabel('c_{m0}');
grid on
axis padded
hold off


figure
plot(p(1,:), alfa_l0(1,:), 'LineWidth', 1.3);
hold on
plot(p(1,:), alfa_l0(2,:), 'LineWidth', 1.3);
plot(p(1,:), alfa_l0(3,:), 'LineWidth', 1.3);
plot(p(1,:), alfa_l0(4,:), 'LineWidth', 1.3);
plot(p(1,:), alfa_l0(5,:), 'LineWidth', 1.3);
plot(p(1,:), alfa_l0(6,:), 'LineWidth', 1.3);
plot(p(1,:), alfa_l0(7,:), 'LineWidth', 1.3);
legend({'f=0.00','f=0.01','f=0.02','f=0.03','f=0.04','f=0.05','f=0.06'}, 'Location','southwest');
xlabel('Posición de camber máxima, p');
ylabel('Ángulo de ataque para sustentación nula, alfa_{l0}');
grid on
axis padded
hold off

