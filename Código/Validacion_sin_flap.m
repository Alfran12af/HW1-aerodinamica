function [Tabla, CL_alfa_grados, alfa_l0] = Validacion_sin_flap(M, f, p, c, xh, eta, U_inf, x_ref)

alfa_limit = 25; % ángulo límite empleado para la zona lineal
alfa_limit = alfa_limit*pi/180; % conversión de grados a rad
alfa_dif = 0.005;
alfa_vector = -alfa_limit:alfa_dif:alfa_limit; % vector alfa

CLDVM = zeros(1, length(alfa_vector));
CM0DVM = zeros(1, length(alfa_vector));

%Matrices CLDVM y CM0DVM para cada ángulo de ataque
for i = 1:length(alfa_vector)
    [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh, eta); %CL y CMLE según DVM
    [G] = Circulacion(M, xcont, xvort, pnorm, alfa_vector(i), U_inf);
    [CLDVM(i), CM0DVM(i)] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa_vector(i), pchord, c, coord);
end

j = 1; % se empieza a leer el vector
fin = 0; % cuando encontremos alfa_l0 se sale del bucle
while (j <= length(alfa_vector) && fin==0)
    if (CLDVM(j) >= 0)
        fin = 1;
    end
    j = j+1;
end

% Cálculo de la pendiente de sustentación
CL_alfa_rad = (CLDVM(j-1)-CLDVM(j))/(alfa_vector(j-1)-alfa_vector(j)); 
CL_alfa_grados = CL_alfa_rad/180*pi; % conversión de 1/rad a 1/grados

% Cálculo del ángulo de sustentación nula
alfa_l0 = alfa_vector(j)-CLDVM(j)/CL_alfa_rad; % interpolación
alfa_l0 = alfa_l0*180/pi; % conversión de rad a grados

% plot(alfa_vector(1,:),CLDVM(1,:)); %plot de CL vs alfa

% Cálculo del coeficiente de momento en el ca
CM0 = CM0DVM(j)-CLDVM(j)*(CM0DVM(j-1)-CM0DVM(j))/(CLDVM(j-1)-CLDVM(j));

% CÁLCULO DEL ERROR RELATIVO:
% Datos tomados experimentalmente - punto 1:
alfa_1 = -6;
alfa_1 = alfa_1*pi/180;
CL_1 = -0.4;
CM0_1 = -0.05;

% Datos tomados experimentalmente - punto 2:
alfa_2 = 2;
alfa_2 = alfa_2*pi/180;
CL_2 = 0.4;
CM0_2 = -0.05;

% Pendiente de sustentación experimental / error relativo:
CL_alfa_exp_rad = (CL_2-CL_1)/(alfa_2-alfa_1);
CL_alfa_exp_grados = CL_alfa_exp_rad/180*pi;
CL_alfa_error_grados = abs(CL_alfa_exp_grados-CL_alfa_grados)/CL_alfa_exp_grados*100;

% Ángulo de sustentación nula experimental / error relativo:
alfa_l0_exp = alfa_1-CL_1/CL_alfa_exp_rad;
alfa_l0_exp = alfa_l0_exp*180/pi;
alfa_l0_error = abs(alfa_l0_exp-alfa_l0)/alfa_l0_exp*100;

% Coeficiente de momento en el ca experimental / error relativo:
CM0_exp = CM0_1-(CM0_2-CM0_1)/(CL_2-CL_1)*CL_1;
CM0_error = abs(CM0_exp-CM0)/CM0_exp*100;

%TABLA RESULTADOS:
Tabla = table({'CL_alfa', 'alfa_l0', 'CM0'}, [CL_alfa_grados alfa_l0 CM0], [CL_alfa_exp_grados alfa_l0_exp CM0_exp], [CL_alfa_error_grados alfa_l0_error CM0_error]);
Tabla.Properties.VariableNames = {'Variable' 'DVM' 'Experimental' 'Error relativo (%)'};





