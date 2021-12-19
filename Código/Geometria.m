function [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh, eta)

%% DISCRETIZACIÓN DE LA GEOMETRÍA
N = M + 1 ;  % Número de puntos en la cuerda
coord = zeros(2,N);   % Matriz vacía para las coordenadas (x,z)
                                 
% Cálculo de las coordenadas de la línea media                              
for i = 1:N
    coord(1,i) = c/2*(1-cos(((i-1)/(N-1))*pi));   % Posición (x) mediante la distribución de cosinus
    % Posición (z) mediante las expresiones de la línea media para perfiles NACA de 4 dígitos
    if ( coord(1,i) < p )  
        coord(2,i) =  c*f/p^2*(2*p*coord(1,i) - coord(1,i).^2);
    else
        coord(2,i) = c*f/(1-p)^2*(1 - 2*p + 2*p*coord(1,i) - coord(1,i).^2); 
    end
end

%% CÁLCULO DEFLEXIÓN FLAP
if ( eta > 0.0 )  % Ángulo de delexión de flap
    % Encontramos el valor de i que hace que x(i) se acerque más a xh:
    diferencia_min = c; % Declaramos inicialmente que la distancia mínima hasta xh sea la longitud de la cuerda
    ih = 1; % Para empezar, la ih es la primera i (1)
    % Hacemos un bucle donde, si la distancia x(i) hasta xh es la más pequeña, almacenamos el valor de i como ih:
    for i = 1 : N 
        diferencia = abs(xh-coord(1,i));
        if (diferencia < diferencia_min)
            diferencia_min = diferencia;
            ih = i;
        end
    end
    
    % Posicón z del hinge en función de la línea media del perfil NACA 4 dígitos
    if ( xh < p ) 
        zh = c*f/p^2*(2*p*coord(1,ih) - coord(1,ih).^2);
    else
        zh = c*f/(1-p)^2*(1 - 2*p + 2*p*coord(1,ih) - coord(1,ih).^2);
    end
    
    Rotacion = [cos(eta) sin(eta); -sin(eta) cos(eta)]; % Matriz de rotación
    coord2 = coord(:, ih:end); % Matriz extraída para x(i) > xh
    coord2(1,:) = coord2(1,:) - coord(1,ih); % Matriz coord2 desplazada al origen (x)
    coord2(2,:) = coord2(2,:) - coord(2,ih); % Matriz coord2 desplazada al origen (z)
    Girados = Rotacion * coord2; % Matriz rotada en el origen
    Girados(1,:) = Girados (1,:) + coord(1,ih); % Matriz rotada en xh
    Girados(2,:) = Girados (2,:) + coord(2,ih); % Matriz rotada en xh
    
    %Reemplazamos los valores calculados en la matriz coord original
    coord(1,ih:end) = Girados(1, :);
    coord(2,ih:end) = Girados(2, :);
end
    
%% CÁLCULO DE LA GEOMETRÏA DE LOS PANELES
pchord = zeros(M,1);
pnorm = zeros(2,M);
ptang = zeros(2,M);
xvort = zeros(2,M);
xcont = zeros(2,M);

for i = 1:M  
    x1 = coord(:,i);  
    x2 = coord(:,i+1);
    
    % Cálculo de la cuerda, el vector normal y el vector tangente de los paneles
    delta_x = x2(1,1) - x1(1,1);
    delta_z = x2(2,1) - x1(2,1);
    pchord(i) = sqrt(delta_x^2+delta_z^2); % Longitud del panel
    pnorm(1,i) = - delta_z / pchord(i) ; % Vectores normales
    pnorm(2,i) = delta_x / pchord(i) ;
    ptang(1,i) = delta_x / pchord(i) ; % Vectores tangenciales
    ptang(2,i) =  delta_z / pchord(i) ; 
    xvort(:,i) = x1(:) + 0.25*pchord(i)*ptang(:,i) ; % Posición del vórtice
    xcont(:,i) = x1(:) + 0.75*pchord(i)*ptang(:,i) ; % Posición del punto de control
end
