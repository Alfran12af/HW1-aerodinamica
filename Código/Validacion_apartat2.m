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

%% CÁLCULO PARA DISTINTAS ALFA (SIN FLAP)

alfa_limit = 30; % ángulo límite empleado para la zona lineal
alfa_dif = 0.005;
alfa = -alfa_limit:alfa_dif:alfa_limit; % vector alfa
alfa = alfa*pi/180; % conversión de grados a rad

Cl_DVM = zeros(1, length(alfa)); % Computed lift coefficient
Cm0_DVM = zeros(1, length(alfa)); % Computed free moment coefficient

% CÁLCULO GEOMETRIA
[coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh, eta);

% Discrete Vortex Method
for i = 1:length(alfa)
    [G, A] = Circulacion(M, xcont, xvort, pnorm, alfa(i), U_inf);
    [Cl_DVM(i), Cm0_DVM(i)] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa(i), pchord, c, coord);
end


%% CÁLCULO ALFA_L0, CL_ALFA (SIN FLAP)

% Cálculo de alfa_l0
i = 0; % Inicio búsqueda
found = 0;

while (i  <= length(alfa)) && (found == 0)
    i = i + 1;
    if Cl_DVM(i) >= 0
        found = 1;
    end
end

x = Cl_DVM(i-1)/(Cl_DVM(i-1)-Cl_DVM(i)); % Interpolació linear
alfa_l0 = x*alfa(i)+(1-x)*alfa(i-1); % alfa_l0

% Cálculo Cl_alfa amb un rang lineal aprox 10º
lin_lim = 10*pi/180; % Linear limit
found = 0; % Binary variable to end search
i = 1; % Initial index for search
i1 = 1; % Index for lower limit of linear range
i2 = 1; % Index for upper limit of linear range

while (i <= length(alfa)) && (found < 2)
    if alfa(i) == -lin_lim
        found = found + 1;
        i1 = i;
    elseif alfa(i) == lin_lim
        found = found + 1;
        i2 = i;
    end
    i = i + 1;
end

Cl_alfa = (Cl_DVM(i2) - Cl_DVM(i1))/(2*lin_lim);



%% CÁLCULO PARA DISTINTOS FLAP-TO CHORD RATIOS

alfa_flap = [0]; % Ángulo de ataque para perfiles con flap
alfa_flap = alfa_flap*pi/180; 

% Cálculo del Cl sin flap
Cl_no_flap = zeros(1, length(alfa_flap));
i = 0;
current = 1; % Current index in alfa_flap

while (i <= length(alfa)-1) && (current <= length(alfa_flap))
    i = i + 1;
    if alfa(i) == alfa_flap(current)
        Cl_no_flap(current) = Cl_DVM(i);
        current = current + 1;
    end
end



% Flap chord ratios, xh y eta
E = 0:0.05:0.4;
xh = c*(1-E);
eta = 10;
eta = eta*pi/180;

% Cálculo Cl con flap

Cl_flap = zeros(length(xh), length(alfa_flap));

for i = 1:length(xh)
    theta_h(i) = pi - acos((xh(i)-0.5)/(c/2));
    alfa_l0_flap(i) = -eta/pi*(pi-theta_h(i)+sin(theta_h(i)));        
    flap_eff_DVM(i, 1) = (alfa_l0_flap(i))/eta;
end



% Flap efficiency corregido
factor = 0.8;
flap_eff_corr = factor*flap_eff_DVM;

% Flap efficiency exp
flap_eff_exp = [0, 0.18, 0.29, 0.38, 0.46, 0.54, 0.60, 0.65, 0.70];


figure
plot(E(1,:), -flap_eff_DVM(:,1));
hold on
plot(E(1,:), -flap_eff_corr(:,1));
plot(E(1,:), flap_eff_exp(1,:));
legend('Eficiencia del flap DVM', 'Eficiencia del flap corregida', 'Eficiencia del flap experimental', 'Location','northwest');
xlabel('Ratio flap-chord, E');
ylabel('Eficiencia del flap');
grid on
hold off