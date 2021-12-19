function [G, A] = Circulacion(M, xcont, xvort, pnorm, alfa, U_inf)

A = zeros(M,M) ;    
RHS = zeros(M,1) ;  

for i = 1:M
    xi = xcont(:,i);  % Posición de los puntos de control
    ni = pnorm(:,i);  % Posición de los vectores normales
    for j = 1:M 
        xv = xvort(:,j); % Posición de los vórtices
        % Cálculo de las velocidades inducidas
        r = sqrt((xi(1,1)-xv(1,1))^2+(xi(2,1)-xv(2,1))^2);
        u = 1/(2*pi)*(xi(2,1)-xv(2,1))/(r^2);
        w = -1/(2*pi)*(xi(1,1)-xv(1,1))/(r^2);
        % Coeficientes de influencia (proyección de la velocidad inducida en ni(:,i))
        A(i,j) = u*ni(1,1)+w*ni(2,1);
    end
    % Proyección de Uinf en ni(:,i)
    RHS(i) = - U_inf * (cos(alfa)*ni(1,1)+sin(alfa)*ni(2,1));
end
G = inv(A)*RHS  ;  % Cálculo de la circulación
