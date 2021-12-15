function [coord, pnorm, ptang, xvort, xcont, pchord] = Geometria(M, f, p, c, xh, eta)

%% DISCRETIZACIÓN DE LA GEOMETRÍA

N = M + 1 ;  % número de puntos en la cuerda

coord = zeros(2,N);   % creación de la matriz vacía para las coordenadas (x,z)
                                                                                                                        
for i = 1:N   % cálculo de las coordenadas de la línea media
    
    coord(1,i) = c/2*(1-cos(((i-1)/(N-1))*pi));   % cálculo de la posición (x) mediante la distribución de cosinus 
    %coord(1,i) = (i-1)/M;
    if ( coord(1,i) < p )  % cálculo de la posición (z) mediante las expresiones de la línea media para perfiles NACA de 4 dígitos
        coord(2,i) =  c*f/p^2*(2*p*coord(1,i) - coord(1,i).^2);
    else
        coord(2,i) = c*f/(1-p)^2*(1 - 2*p + 2*p*coord(1,i) - coord(1,i).^2); 
    end
     
end


%% CÁLCULO DEFLEXIÓN FLAP:

if ( eta > 0.0 )  %Angulo de delexin de Flap
    
    % Encontramos el valor de i que hace que x(i) se acerque más a xh:
    diferencia_min = c; % Declaramos inicialmente que la distancia mínima hasta xh sea la longitud de la cuerda
    ih = 1; % Para empezar, la ih es la primera i (1)
    
    %Hacemos un bucle donde, si la distancia x(i) hasta xh es la más
    %pequeña, almacenamos el valor de i como ih:
    for i = 1 : N 
        diferencia = abs(xh-coord(1,i));
        if (diferencia < diferencia_min)
            diferencia_min = diferencia;
            ih = i;
        end
    end
            
    if ( xh < p )  % hinge z position according to NACA4 mean line function
        zh = c*f/p^2*(2*p*coord(1,ih) - coord(1,ih).^2);
    else
        zh = c*f/(1-p)^2*(1 - 2*p + 2*p*coord(1,ih) - coord(1,ih).^2);
    end
    
    Rotacion = [cos(eta) sin(eta); -sin(eta) cos(eta)]; % Matriz de rotación
    coord2 = coord(:, ih:end); % Matriz extraída para x(i) > xh
    coord2(1,:) = coord2(1,:) - coord(1,ih); % Matriz coord2 desplazada al origen (x)
    coord2(2,:) = coord2(2,:) - coord(2,ih); % Matriz coord2 desplazada al origen (z)
   
    Girados = Rotacion * coord2; % Matriz en el origen rotada
    
    Girados(1,:) = Girados (1,:) + coord(1,ih); % Matriz rotada en xh
    Girados(2,:) = Girados (2,:) + coord(2,ih); % Matriz rotada en xh
    
    %Reemplazamos los valores calculados en la matriz coord original
    coord(1,ih:end) = Girados(1, :);
    coord(2,ih:end) = Girados(2, :);
    
end
    
            
       

pchord = zeros(M,1) ;   % allocate storage space
pnorm = zeros(2,M) ;
ptang = zeros(2,M) ;
xvort = zeros(2,M) ;
xcont = zeros(2,M) ;

% CÁLCULO DE LA GEOMETRÍA DE LOS PANELES

for i = 1:M  
    x1 = coord(:,i);  
    x2 = coord(:,i+1);
    
    % calculate panel's chord, normal and tangent vectors
    
    delta_x = x2(1,1) - x1(1,1);
    delta_z = x2(2,1) - x1(2,1);
    pchord(i) = sqrt(delta_x^2+delta_z^2);
    pnorm(1,i) = - delta_z / pchord(i) ; % components of the normal vector
    pnorm(2,i) = delta_x / pchord(i) ;
    ptang(1,i) = delta_x / pchord(i) ; % components of the tangent vector
    ptang(2,i) =  delta_z / pchord(i) ; 
    
    % calculate position of the vortex and control point
    
    xvort(:,i) = x1(:) + 0.25*pchord(i)*ptang(:,i) ;
    xcont(:,i) = x1(:) + 0.75*pchord(i)*ptang(:,i) ;
    
end

%{
% Comprovación gráfica de los resultados
figure(1)

plot(coord(1,:),coord(2,:), '-o');
hold on
quiver(xcont(1,:), xcont(2,:), pnorm(1,:), pnorm(2,:));
axis padded
hold off
%}


