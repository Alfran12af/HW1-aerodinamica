clear all
clc
% CODI FINAL!!!

%% DATOS
M =  200;  % número de paneles
f =  0;  % curvatura máxima 
p =  0;  % posición curvatura máxima
c = 1; %longitud cuerda
U_inf = 100; % velocidad de corriente libre
x_ref = c/4; % referencia: c/4



%% Variacion de Cl y Cm para distintas posiciones del Flap alfa=0
% Flap chord ratios, xh y eta
E =  [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4];
xh = c*(1-E);
eta = linspace(0, 45, 20);
eta = eta*pi/180;

% Cálculo Cl con flap
EF = zeros (size(xh));
Clf_pediente= zeros(1,length(xh));
CMf_pendiente = zeros(1,length(xh));


for j = 1:length(xh)
    IncAlpha0 = zeros(1,length(eta));
    Cl0_AUX = zeros(1,length(eta));
    Clf_Vec_AUX = zeros(1,length(eta));
    CMf_Vec_AUX = zeros(1,length(eta));
    
    for i = 1:length(eta)
        % Cálculo de los coeficientes según el método DVM para cada eta
        [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh(j), eta(i));
        [G, A] = Circulacion(M, xcont, xvort, pnorm, 0, U_inf);
        [Cl_DVM_AUX, CM_DVM_AUX,Clf_DVM_AUX,CMf_DVM_AUX] = CoeficientesDVM(M, U_inf, G, xvort, xh(j), 0, pchord, c, coord,xh(j));
        Clf_Vec_AUX(i) =  Clf_DVM_AUX; 
        CMf_Vec_AUX(i) = CMf_DVM_AUX;
        
        disp(['eta ' num2str(Clf_DVM_AUX)]);
    end
    % Regresiones lineales para obtener cada las pendientes de cada coeficientes
    line1 = polyfit(eta,Clf_Vec_AUX,1);
    line2 = polyfit(eta,CMf_Vec_AUX,1);
    Clf_pediente(j) = line1(1);
    CMf_pendiente(j) = line2(1);
    disp(['EF' num2str(j)]);
end


q = 1/2*rho*(U_inf)^2;

% Cálculo de los vectores de Clf, CM, L, M
for i = 1:length(xh)
    for j = 1:length(eta_v)
            delta_CL (i,j) = Clf_pediente(i)*eta_v(j);
            delta_CM0 (i, j) = CMf_pediente(i)*eta_v(j);
            L (i, j) = q*delta_CL(i,j)*(1-xh(i))*b;
            M0 (i, j) = q*delta_CM0(i,j)*(1-xh(i))^2*b;
    end
end

%% FIGURAS

% Declaración de rectas para conseguir la leyenda buscada
L(8,:) = L(7,:)*1.01;
L(9,:) = L(7,:)*1.01;
M0(8,:) = M0(7,:)*1.01;
M0(9,:) = M0(7,:)*1.01;
delta_CL(8,:) = delta_CL(7,:)*1.01;
delta_CL(9,:) = delta_CL(7,:)*1.01;
delta_CM0(8,:) = delta_CM0(7,:)*1.01;
delta_CM0(9,:) = delta_CM0(7,:)*1.01;

% Declaración de colores para todas las variables
str1 = '#A2142F';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
str2 = '#77AC30';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
str3 = '#D95319';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;
str4 = '#77AC30';
color4 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;


figure
yyaxis left
plot(eta_v(1,:), delta_CL(4,:), '-b');
ylabel('Variación del coeficiente de sustentación, C_{l_{flap}}');
hold on
plot(eta_v(1,:), delta_CL(5,:), '-g');
plot(eta_v(1,:), delta_CL(6,:), '-', 'Color', color3);
plot(eta_v(1,:), delta_CL(7,:), '-r');
plot(eta_v(1,:), delta_CL(8,:), '-w');
plot(eta_v(1,:), delta_CL(9,:), '-w');
axis padded

yyaxis right

ylabel('Variación del coeficiente de momento de cabeceo, C_{m_{hinge}}');


plot(eta_v(1,:), delta_CM0(4,:), '-.b');
hold on
plot(eta_v(1,:), delta_CM0(5,:), '-.g');
plot(eta_v(1,:), delta_CM0(6,:), '-.', 'Color', color3);
plot(eta_v(1,:), delta_CM0(7,:), '-.r');
title('Análisis de \Delta C_{l_{flap}} y \Delta C_{m_{hinge}} en función de \eta y xh/c');
legend({'xh/c = 0.85','xh/c = 0.80','xh/c = 0.75','xh/c = 0.70', '\Delta C_{l_{flap}}', '\Delta C_{m_{hinge}}'}, 'Location','west');
xlabel('Ángulo de deflexión, \eta [rad]');
grid on
axis padded
hold off


figure
yyaxis left
plot(eta_v(1,:), L(4,:), '-b');
ylabel('Variación de la fuerza de sustentación, L_{flap} [N/m]');
hold on
plot(eta_v(1,:), L(5,:), '-g');
plot(eta_v(1,:), L(6,:), '-','Color', color1);
plot(eta_v(1,:), L(7,:), '-r');
plot(eta_v(1,:), L(8,:), '-w');
plot(eta_v(1,:), L(9,:), '-w');
axis padded

yyaxis right

ylabel('Variación del momento de cabeceo, M_{hinge} [N·m/m]');

plot(eta_v(1,:), M0(4,:), '-.b');
hold on
plot(eta_v(1,:), M0(5,:), '-.g');
plot(eta_v(1,:), M0(6,:), '-.','Color', color1);
plot(eta_v(1,:), M0(7,:), '-.r');

title('Análisis de L_{flap} y M_{hinge} por unidad de envergadura en función de \eta y xh/c');
legend({'xh/c = 0.85','xh/c = 0.80','xh/c = 0.75','xh/c = 0.70', '\Delta L_{flap}', '\Delta M_{hinge}'}, 'Location','west');
xlabel('Ángulo de deflexión, \eta [rad]');
grid on
axis padded
hold off



