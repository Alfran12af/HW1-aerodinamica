function [G, A] = Circulacion(M, xcont, xvort, pnorm, alfa, U_inf)


%% (2) Assembly and solution of the equations system
% 


A = zeros(M,M) ;  % allocate variables
RHS = zeros(M,1) ;  

for i = 1:M  % scan control points
    
    xi = xcont(:,i)   ;  % control point's position
    ni = pnorm(:,i)   ;  % normal vector
    
    for j = 1:M  % scan vortices
        
        xv = xvort(:,j)  ;  % vortex position
        r = sqrt((xi(1,1)-xv(1,1))^2+(xi(2,1)-xv(2,1))^2);
        u = 1/(2*pi)*(xi(2,1)-xv(2,1))/(r^2)   ; % induced velocities
        w = -1/(2*pi)*(xi(1,1)-xv(1,1))/(r^2) ;
            
        A(i,j) = u*ni(1,1)+w*ni(2,1);  % influence coefficient (projection of the induced velocity onto ni(:,i))
            
    end
    
    RHS(i) = - U_inf * (cos(alfa)*ni(1,1)+sin(alfa)*ni(2,1))  ;  % projection of Uinf onto ni(:,i), don't forget the minus sign!
        
end

G = inv(A)*RHS  ;  % circulation
