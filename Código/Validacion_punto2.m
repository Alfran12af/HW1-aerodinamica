clear all
clc
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

%% CÁLCULO PARA DISTINTAS ALFA (SIN FLAP)

alfa_limit = 10; % ángulo límite empleado para la zona lineal
alfa_dif = 0.01;
alfa = -alfa_limit:alfa_dif:alfa_limit; % vector alfa
alfa = alfa*pi/180; % conversión de grados a rad

Cl_DVM = zeros(1, length(alfa)); % cálculo del coeficiente de sustentación por DVM
Cm0_DVM = zeros(1, length(alfa)); % cálculo del coeficiente de momento libre por DVM

% CÁLCULO GEOMETRIA
[coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh, eta);

% Discrete Vortex Method
for i = 1:length(alfa)
    [G, A] = Circulacion(M, xcont, xvort, pnorm, alfa(i), U_inf);
    [Cl_DVM(i), Cm0_DVM(i)] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa(i), pchord, c, coord,xh);
end

%% CÁLCULO ALFA_L0, CL_ALFA (SIN FLAP)

% Cálculo de alfa_l0
i = 0; % Inicio búsqueda
found = 0;

% Búsqueda de la sustentación nula
while (i  <= length(alfa)) && (found == 0)
    i = i + 1;
    if Cl_DVM(i) >= 0
        found = 1;
    end
end

x = Cl_DVM(i-1)/(Cl_DVM(i-1)-Cl_DVM(i)); % Interpolación linear
alfa_l0 = x*alfa(i)+(1-x)*alfa(i-1); % Cálculo del ángulo de sustentación nula

% Cálculo Cl_alfa en el rango lineal 8º
lin_lim = 8*pi/180; % Ángulo del límite lineal
found = 0;
i = 1; % Inicio de búsqueda
i1 = 1; % Límite inferior del rango lineal
i2 = 1; % Límite superior del rango lineal

while (i <= length(alfa)) && (found < 2)
    if alfa(i) == -lin_lim % Búsqueda del límite inferior
        found = found + 1;
        i1 = i;
    elseif alfa(i) == lin_lim % Búsqueda del límite superior
        found = found + 1;
        i2 = i;
    end
    i = i + 1;
end

Cl_alfa = (Cl_DVM(i2) - Cl_DVM(i1))/(2*lin_lim); % Cálculo de la pendiente de la curva de sustentación



%% CÁLCULO PARA DISTINTOS FLAP-TO CHORD RATIOS

alfa_flap = [0]; % Ángulo de ataque para perfiles con flap
alfa_flap = alfa_flap*pi/180; 

% Cálculo del Cl sin flap
Cl_no_flap = zeros(1, length(alfa_flap));
i = 0;
current = 1;

while (i <= length(alfa)-1) && (current <= length(alfa_flap))
    i = i + 1;
    if alfa(i) == alfa_flap(current)
        Cl_no_flap(current) = Cl_DVM(i);
        current = current + 1;
    end
end


% Flap chord ratios, xh y eta
E =  [0.15 0.2 0.25 0.3];
xh = c*(1-E);
eta = linspace(0, 45, 10);
eta = eta*pi/180;

% Cálculo Cl con flap
EF = zeros (size(xh));

for j = 1:length(xh)
    % Declaración de los vectores de cada variable
    IncAlpha0 = zeros(1,length(eta));
    Cl0_AUX = zeros(1,length(eta));
    for i = 1:length(eta)
        % Cálculo de la geometría para cada eta
        [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh(j), eta(i));
        Cl_DVM_AUX = zeros(1,length(alfa));
        
        for k = 1:length(alfa)
            % Cálculo de los coeficientes DVM para cada alfa
            [G, A] = Circulacion(M, xcont, xvort, pnorm, alfa(k), U_inf);
            [Cl_DVM_AUX(k), t1,t2,t3] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa(k), pchord, c, coord,1);
        end
        % Regresión lineal para obtener el ángulo de sustentación nula
        line = polyfit(alfa,Cl_DVM_AUX, 1); 
        alfa0 = -line(2)/line(1);
        Cl0_AUX(i) = alfa0;
        disp(['eta ' num2str(i)]);
    end
    % Regresión lineal para obtener el factor de eficiencia del flap
    line = polyfit(eta,Cl0_AUX,1);
    EF(j) = line(1); 
    disp(['EF' num2str(j)]);
end


% Valores del factor de eficiencia experimental
flap_eff_exp = [0.38, 0.46, 0.54, 0.60];

% Vector con posibles factores de corrección 
factor = 0.8:0.001:0.9;
largo = length(factor);
dif = zeros(1,largo);

% Cálculo del factor de eficiencia del flap corregido

for i = 1:largo % Algoritmo para calcular el factor de correccion con menos error
   for j = 4:7
      dif(i) = dif(i) + abs(flap_eff_exp(j)+factor(i)*EF(j,1));
   end
end

dif_min = min(dif);
num = find(dif==dif_min);
factor_2 = factor(num);

flap_eff_corr = factor_2*EF; % Cálculo del factor del flap corregido con el factor de corrección óptimo

%% FIGURAS

figure(1)
plot(E, -EF), 'r';
hold on
plot(E, -flap_eff_corr,'g');
plot(E, flap_eff_exp,'b');
legend('Eficiencia del flap DVM', 'Eficiencia del flap corregida', 'Eficiencia del flap experimental', 'Location','northwest');
xlabel('Ratio flap-chord, E');
ylabel('Eficiencia del flap');
grid on
hold off