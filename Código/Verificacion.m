clear all

M =  200;  % número de paneles
f =  0.02;  % curvatura máxima 
p =  0.4;  % posición curvatura máxima
c = 1; %longitud cuerda

alfa = 4;% ángulo de ataque del análisis
alfa = alfa*pi/180; % conversión de grados a rad

xh =  1.0;  % posición hinge 
eta =  0.0; % ángulo deflexión flap

U_inf = 1; % velocidad de corriente libre
x_ref = 0; % referencia: LE

Matriu_CL = zeros(4,M);
Matriu_CMLE = zeros(4,M);

for i=1:M
    Matriu_CL(1, i) = i;
    Matriu_CMLE(1, i) = i;
    
    % Cálculo de los coeficientes según el método TAT
    [CLTAT, CMLETAT] = CoeficientesTAT(p, f, c, alfa); %CL y CMLE según TAT
    Matriu_CL(2, i) = CLTAT;
    Matriu_CMLE(2, i) = CMLETAT;
    
    % Cálculo de los coeficientes según el método DVM
    [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(i, f, p, c, xh, eta);
    [G] = Circulacion(i, xcont, xvort, pnorm, alfa, U_inf);
    [CLDVM, CMLEDVM] = CoeficientesDVM(i, U_inf, G, xvort, x_ref, alfa, pchord, c, coord);
    Matriu_CL(3, i) = CLDVM;
    Matriu_CMLE(3, i) = CMLEDVM;
    
    % Cálculo del error relativo
    Matriu_CL(4,i) = abs(abs(Matriu_CL(3, i)-Matriu_CL(2, i))/Matriu_CL(2, i))*100;
    Matriu_CMLE(4,i) = -abs(abs(Matriu_CMLE(3, i)-Matriu_CMLE(2, i))/Matriu_CMLE(2, i))*100;
end

%% FIGURAS

figure
yyaxis left
plot(Matriu_CL(1,:), Matriu_CL(2,:), 'b','LineWidth',1.3);
xlabel('Número de paneles');
ylabel('Coeficiente de sustentación C_L');
title('Coeficiente de sustentación C_L según TAT y DVM')
hold on
plot(Matriu_CL(1,:), Matriu_CL(3,:), 'g -','LineWidth',1.3);
yyaxis right
ylabel('Error relativo C_L (%)')
plot(Matriu_CL(1,:), Matriu_CL(4,:), 'r','LineWidth',1.3);
grid on
ylim([-3 35])
legend('C_L según TAT', 'C_L según DMV', 'Error relativo C_L (%)', 'Location', 'East');
hold off


figure
yyaxis left
plot(Matriu_CMLE(1,:), Matriu_CMLE(2,:), 'b','LineWidth',1.3);
xlabel('Número de paneles');
ylabel('Coeficiente de momento C_{MLE}');
ylim([-0.27 0])
title('Coeficiente de momento C_{MLE} según TAT y DVM')
hold on
plot(Matriu_CMLE(1,:), Matriu_CMLE(3,:), 'g -','LineWidth',1.3);
yyaxis right
ylabel('Error relativo C_{MLE} (%)')
plot(Matriu_CMLE(1,:), Matriu_CMLE(4,:), 'r','LineWidth',1.3);
ylim([-55 3])
grid on
legend('C_{MLE} según TAT', 'C_{MLE} según DMV', 'Error relativo C_{MLE} (%)', 'Location', 'East');
hold off

figure
plot(coord(1,:),coord(2,:), '-o');
hold on
quiver(xcont(1,:), xcont(2,:), pnorm(1,:), pnorm(2,:));
quiver(xcont(1,:), xcont(2,:), ptang(1,:), ptang(2,:));
axis padded
hold off
