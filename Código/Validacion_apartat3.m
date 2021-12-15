clear all
clc
% CODI FINAL!!!

%% DATOS
M =  200;  % número de paneles
f =  0.0;  % curvatura máxima 
p =  0.0;  % posición curvatura máxima
c = 1; %longitud cuerda
U_inf = 100; % velocidad de corriente libre
x_ref = c/4; % referencia: c/4
rho = 1.225;

E = 0:0.05:0.4;
xh = c*(1-E);

eta_v = 0:1:20;

eta_v = eta_v*pi/180;

%% VARIACION FLAP CL Y CM0; FORCE AND PITCHING MOMENT ACTING ON HINGE

q = 1/2*rho*(U_inf)^2;

for i = 1:length(xh)
    for j = 1:length(eta_v)
            theta_h(i, j) = pi - acos((xh(i)-0.5)/(c/2));
    
            delta_CL (i, j) = 2*(pi-theta_h(i)+sin(theta_h(i)))*eta_v(j);
            delta_CM0 (i, j) = -1/2*sin(theta_h(i))*(1-cos(theta_h(i)))*eta_v(j);
            L (i, j) = q*delta_CL(i,j)*xh(i);
            M0 (i, j) = q*delta_CM0(i,j)*xh(i)*c;
    end
end


str1 = '#A2142F';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
str2 = '#77AC30';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
str3 = '#D95319';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;



figure
yyaxis left
plot(eta_v(1,:), delta_CL(1,:), ':+r');
ylabel('Variación del coeficiente de sustentación, C_{L}');
hold on
plot(eta_v(1,:), delta_CL(2,:), ':+g');
plot(eta_v(1,:), delta_CL(3,:), ':+b');
plot(eta_v(1,:), delta_CL(4,:), ':+c');
plot(eta_v(1,:), delta_CL(5,:), ':+m');
plot(eta_v(1,:), delta_CL(6,:), ':+y');
plot(eta_v(1,:), delta_CL(7,:), ':+','Color', color1);
plot(eta_v(1,:), delta_CL(8,:), ':+','Color', color2);
plot(eta_v(1,:), delta_CL(9,:), ':+', 'Color', color3);

axis padded

yyaxis right

ylabel('Variación del coeficiente de momento de cabeceo, C_{M0}');


plot(eta_v(1,:), delta_CM0(1,:), '-r');
hold on
plot(eta_v(1,:), delta_CM0(2,:), '-g');
plot(eta_v(1,:), delta_CM0(3,:), '-b');
plot(eta_v(1,:), delta_CM0(4,:), '-c');
plot(eta_v(1,:), delta_CM0(5,:), '-m');
plot(eta_v(1,:), delta_CM0(6,:), '-y');
plot(eta_v(1,:), delta_CM0(7,:), '-', 'Color', color1);
plot(eta_v(1,:), delta_CM0(8,:), '-', 'Color', color2);
plot(eta_v(1,:), delta_CM0(9,:), '-', 'Color', color3);

legend({'xh = 1','xh = 2','xh = 3','xh = 4','xh = 5','xh = 6','xh = 7','xh = 8','xh = 9'}, 'Location','northwest');
xlabel('Ángulo de deflacción, eta');
plot(0.3, 0, 'w', 'DisplayName','- \Delta C_{M0}');
plot(0.3, 0, 'w','DisplayName','+ \Delta C_L');
grid on
axis padded
hold off






figure
yyaxis left
plot(eta_v(1,:), L(1,:), ':+r');
ylabel('Variación de la fuerza de sustentación, L');
hold on
plot(eta_v(1,:), L(2,:), ':+g');
plot(eta_v(1,:), L(3,:), ':+b');
plot(eta_v(1,:), L(4,:), ':+c');
plot(eta_v(1,:), L(5,:), ':+m');
plot(eta_v(1,:), L(6,:), ':+y');
plot(eta_v(1,:), L(7,:), ':+','Color', color1);
plot(eta_v(1,:), L(8,:), ':+','Color', color2);
plot(eta_v(1,:), L(9,:), ':+', 'Color', color3);

axis padded

yyaxis right

ylabel('Variación del momento de cabeceo, M_0');


plot(eta_v(1,:), M0(1,:), '-r');
hold on
plot(eta_v(1,:), M0(2,:), '-g');
plot(eta_v(1,:), M0(3,:), '-b');
plot(eta_v(1,:), M0(4,:), '-c');
plot(eta_v(1,:), M0(5,:), '-m');
plot(eta_v(1,:), M0(6,:), '-y');
plot(eta_v(1,:), M0(7,:), '-', 'Color', color1);
plot(eta_v(1,:), M0(8,:), '-', 'Color', color2);
plot(eta_v(1,:), M0(9,:), '-', 'Color', color3);

legend({'xh = 1','xh = 2','xh = 3','xh = 4','xh = 5','xh = 6','xh = 7','xh = 8','xh = 9'},'Location','northwest');
xlabel('Ángulo de deflacción, eta');
plot(0.3, 0, 'w', 'DisplayName','- \Delta M_0');
plot(0.3, 0, 'w','DisplayName','+ \Delta  L');

grid on
axis padded
hold off




