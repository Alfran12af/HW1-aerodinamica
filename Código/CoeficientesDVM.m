function [CL, CMLE] = CoeficientesDVM(M, U_inf, G, xvort, x_ref, alfa, pchord, c, coord)


%% (3) Loads calculation
% 

CL = 0.0 ; CMLE = 0.0  ;  % set to zero
VectorCp = zeros(1,M);
for i = 1:M   % scan panels' vortices
    
    delta_CL = 2/(U_inf*c)*G(i,1)   ; %panel's lift using K-J
    delta_CM = -2/(U_inf*c^2)*G(i,1)*(xvort(1,i)-x_ref)*cos(alfa)  ; %moment, definim
    delta_CP = (2*G(i,1))/(U_inf*pchord(i,1));
    % you can also compute delta_CP in this loop, in this case
    % allocate memory space in advance, i.e deltaCP = zeros(M,1)
        
    CL = CL + delta_CL ;   % accumulate panels' contributions
    CMLE = CMLE + delta_CM  ;
    VectorCp(1,i)= delta_CP;
    
end

X_med = zeros (M:1);
for i=1:M
    X_med(i,1)=(coord(1,i+1)+coord(1,i))/2;
end


% You can plot the distribution of circulation and deltaCP for further analysis
% Flap lift and moment about the hinge can be similarly obtained but using only the vortex along the flap and taking moments about the hinge position (xh,zh)

% plot of the distribution of circulation and deltaCP

%{
figure
hold on
plot(X_med, VectorCp(1,:));
grid on
grid minor
hold off

figure
hold on
bar(X_med, G(:,1));
plot(X_med, G(:,1),'g', 'LineWidth', 2.5);
grid on
grid minor
hold off
%}
