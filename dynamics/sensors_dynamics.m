%% 

% Program to solve the Linblad master equation for COUPLING OF A TWO LEVEL
% ATOM WITH A MODE

tic

% INITIAL STATE

% Initial state of the atom = EXCITED STATE with some COHERENCE

rho_at_in = [0 -0.5;0.5 1];

% Initial state of the mode = GROUND STATE

n = 3;
rho_mode_in = zeros(n);
rho_mode_in(1,1) = 1;

% Initial state of the sensor = GROUND STATE

rho_sens_in = [1 -0.5;0.5 0];

% total initial rho

rho_bath = kron(rho_at_in,rho_mode_in);

rhotot_in = kron(rho_sens_in,rho_bath);

rhotot_in_vector = reshape(rhotot_in,[],1);

%% Solution of master equation

[t,y] = ode45(@func_sensors,[0 100],rhotot_in_vector);

    % the initial condition are derived by the expression of rho_th
    % in the main function
   
y2 = reshape(y,[size(y,1),sqrt(size(y,2)),sqrt(size(y,2))]);
    % it's a solid figure tx6x6

y3 = zeros(size(y,1),2,2);

for time = 1:size(y2,1)
   
    temporary = squeeze(y2(time,:,:));     
        %temporary is a matrix 6x6, calculated at each istant of time
        
    y3(time,:,:) = TrX(temporary,2,[length(rho_sens_in),length(rho_bath)]); 
    
        %y3 is a solid figure (timex2x2), calculated at each istant of time
    
end

%%
figure; 

plot(t,y3(:,1,1),'-',t,y3(:,1,2),'-',t,y3(:,2,1),'-',t,y3(:,2,2),'-')

legend('\rho_{11}', '\rho_{12}','\rho_{21}','\rho_{22}',...
    'Location','bestoutside');

title('Dynamics of the sensor',...
    'color','k','fontsize', 18,'fontname','helvetica',...
    'fontunits','normalized','fontweight','normal');

xlabel('t(ps)','color','k','fontsize', 12);
ylabel('\rho(t)','color','k','fontsize', 12);

toc