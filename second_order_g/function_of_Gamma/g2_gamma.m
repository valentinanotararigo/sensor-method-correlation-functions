
tic

% operators

n=7;   % to calculate sigma_z I need a_tmp
ntrunc = [1,n-1,1,1]; 
    %the two level atom, the mode, 2 senosors

sigma_m = func_operators(ntrunc,1);
a = func_operators(ntrunc,2);
s1 = func_operators(ntrunc,3);
s2 = func_operators(ntrunc,4);

sigma_p = sigma_m';

% parameters

w_mode = 0;                
ht=1;
g = 1;

gamma_atom = 0.01*g;
gamma_pump = gamma_atom;
gamma_mode = 0.1*g;

epsilon1 = 0.00001;
epsilon2 = 0.00001;

R = sqrt(g^2-((gamma_mode-gamma_atom)/4)^2);
R2_minus = sqrt(2*g^2-((gamma_mode-gamma_atom)/4)^2)...
           -sqrt(g^2-((gamma_mode-gamma_atom)/4)^2);

w_2 = R;
w_sens = [R2_minus; -R];

nloop = 100;                                                                           ;
gamma_sens = logspace(-4,2,nloop);
% -4 and 2 are the powers of 10, so in a linear scale, 
% I would have: linspace(10^(-4),100,nloop)
% An alternative way could be:
% tmp = linspace(log(10^(-4)),log(100),nloop);
% gamma_sens = exp(tmp);  

%% steady state

H_JC = ht * (sigma_p*a + sigma_m*a');
H_int = ht * (epsilon1 * (a*s1'+a'*s1) + epsilon2 * (a*s2'+a'*s2));
H_inv = H_JC + H_int;

L_atom = dissipator(gamma_atom, sigma_m, H_inv);
L_pump = dissipator(gamma_pump, sigma_p, H_inv);  
L_mode = dissipator(gamma_mode, a, H_inv);

L_inv = L_atom + L_pump + L_mode;

for j = 1:length(w_sens)
    
    H_sens = ht * (w_sens(j) * s1' * s1 + w_2 * s2' * s2);
    H = H_inv + H_sens;
    
    LH = -1i * (kron(H.',eye(length(H)))-kron(eye(length(H)),H));

    for k = 1:nloop
          
        L_sens = dissipator(gamma_sens(k), s1, H_inv) +...
            dissipator(gamma_sens(k), s2, H_inv); 
        L = LH + L_inv + L_sens;
     
        if k==1
            [A,D] = eigs(L,1,0);    
            %eigs(A,k,sigma) returns k eigenvalues with value=sigma
        else
            [A,D] = eigs(L,1,0,opts);
            %eigs(A,K,sigma,opts) specifies an options structure. 
        end
    
        rho_ss_tmp(:,j,k) = A;
        opts.v0 = A;           %opts.v0=starting vector

    end
end

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)),size(rho_ss_tmp,2),size(rho_ss_tmp,3)]);

%% g2(0) as a function of the sensor linewidth Gamma

rho_ss_norm = zeros(size(rho_ss,1),size(rho_ss,2),...
                length(w_sens),nloop);
rho_corr1_tmp = zeros(size(rho_ss,1),size(rho_ss,2),...
                length(w_sens),nloop);
rho_corr1 = zeros(length(w_sens),nloop);
n1 = zeros(length(w_sens),nloop);
rho_corr2_tmp = zeros(size(rho_ss,1),size(rho_ss,2),...
                length(w_sens),nloop);
rho_corr2 = zeros(length(w_sens),nloop);
n2 = zeros(length(w_sens),nloop);
rho_2ndorder_tmp = zeros(size(rho_ss,1),size(rho_ss,2),...
                   length(w_sens),nloop);
rho_2ndorder = zeros(length(w_sens),nloop);
g2 = zeros(length(w_sens),nloop);
g2_norm = zeros(length(w_sens),nloop);

for j = 1:length(w_sens)
    
    for k = 1:nloop
        
        rho_ss_norm(:,:,j,k) = rho_ss(:,:,j,k) / trace(rho_ss(:,:,j,k));
        
        rho_corr1_tmp(:,:,j,k) = s1' * s1 * rho_ss_norm(:,:,j,k);
        rho_corr1(j,k) = trace(rho_corr1_tmp(:,:,j,k));
        n1(j,k) = gamma_sens(k)/(2*pi*epsilon1^2) * rho_corr1(j,k);
        
        rho_corr2_tmp(:,:,j,k) = s2' * s2 * rho_ss_norm(:,:,j,k);
        rho_corr2(j,k) = trace(rho_corr2_tmp(:,:,j,k));
        n2(j,k) = gamma_sens(k)/(2*pi*epsilon2^2) * rho_corr2(j,k);
    
        rho_2ndorder_tmp(:,:,j,k) = s1' * s1 * s2' * s2 * rho_ss_norm(:,:,j,k);
        rho_2ndorder(j,k) = trace(rho_2ndorder_tmp(:,:,j,k));
        
        g2(j,k) = (gamma_sens(k)*gamma_sens(k))/...
            (2*pi*epsilon1*epsilon2)^2 * rho_2ndorder(j,k);
        g2_norm(j,k) = g2(j,k) / (n1(j,k)*n2(j,k));
    end
end

%% figure

u = [0.055 0.21 0.41 1.9995];
%u = [0.06 0.4 0.55 1.9995];

figure;   %then modify in loglog scale
hold on;

colorspec = {[0 0.5 0]; [1 0.8 0]};  % green,yellow

for j = 1:length(w_sens)
    
    line(gamma_sens, g2_norm(j,:),'Color', colorspec{j});
    adiff = gca; % current axes
    adiff_pos = adiff.Position; % position of first axes
    axis([10^(-4) 100 0 6]);
    
    %plot(gamma_sens,g2_norm(j,:));
    
    legend('\omega_1=R_2^-','\omega_1=-R',...
        'fontsize', 12,'Location','best');
end

xlabel('(\Gamma/g','color','k','fontsize', 20);
ylabel('g^{(2)}(\Gamma)','color','k','fontsize', 20);

hold on;

au = axes('Position',adiff_pos,...
            'XAxisLocation','top','YAxisLocation','right',...
            'XTick',u,'Color','none');
%line(u,ones(1,length(u)),'Color','none','Parent',au);
axis([10^(-4) 100 0 6]);

toc
    
    
    
    
    
    
    
