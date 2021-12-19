function [CL, CMLE] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa, pchord, c, coord)

% Iniciamos con los coeficientes a cero
CL = 0.0; CMLE = 0.0;
VectorCp = zeros(1,M);

for i = 1:M   % Análisis panel a panel
    delta_CL = 2/(U_inf*c)*G(i,1); % Cálculo del coeficiente de sustentación por Kutta-Joukowski
    delta_CM = -2/(U_inf*c^2)*G(i,1)*(xvort(1,i)-x_ref)*cos(alfa); % Cálculo del coeficiente de momento en el borde de ataque
    delta_CP = (2*G(i,1))/(U_inf*pchord(i,1)); % Cálculo del coeficiente de presión
    
    % Suma de las contribuciones de cada panel
    CL = CL + delta_CL;
    CMLE = CMLE + delta_CM;
    VectorCp(1,i)= delta_CP;
end

X_med = zeros (M:1);
for i=1:M
   % Punto x medio de cada panel
    X_med(i,1)=(coord(1,i+1)+coord(1,i))/2;
end
