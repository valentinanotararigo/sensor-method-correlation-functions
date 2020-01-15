%% Calculation of the 2 steady states for configurations ii and iv

tic

clear;

rho_ss_tmp = func_ss;

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)), size(rho_ss_tmp,2)]);

load parameters.mat;
       
n_configurations = 2;
    
%% Calculation of <n1(0)>

load my_operators.mat n ntrunc s1 s2;

    % I'm loading only some of the operators stored in the file
    % my_operators.mat

for j = 1:n_configurations
    
n1_mean_tmp(:,:,j) = s1' * s1 * rho_ss(:,:,j)/trace(rho_ss(:,:,j));
n1_mean(j) = trace(n1_mean_tmp(:,:,j));
S_n1(j) = gamma_sens/(2*pi*epsilon1^2) * n1_mean(j);

end

%% Dynamics of n2(t) with initial condition n2(0)

% w1=R2^-

n2 = s2'* s2;
n2_in_vector = reshape(n2,[],1);

[T1,n2_T1] = ode45(@func_n2_T1,[0 60],n2_in_vector);
%[T1,n2_T1] = ode45(@func_n2_T1,[-60 0],n2_in_vector);

n2_bis = reshape(n2_T1,[size(n2_T1,1),sqrt(size(n2_T1,2)),...
    sqrt(size(n2_T1,2))]);    % it's a solid figure tx56x56
     
for time = 1:size(n2_bis,1)
   
    n2_f1 = permute(n2_bis,[2 3 1]);  % it's a solid figure 56x56xt
end

% w1=R

[T2,n2_T2] = ode45(@func_n2_T2,[0 60],n2_in_vector);
%[T2,n2_T2] = ode45(@func_n2_T2,[-60 0],n2_in_vector);

n2_T2_bis = reshape(n2_T2,[size(n2_T2,1),sqrt(size(n2_T2,2)),...
    sqrt(size(n2_T2,2))]);
     
for time2 = 1:size(n2_T2_bis,1)
   
    n2_f2 = permute(n2_T2_bis,[2 3 1]);
end

%% Calculation of <n2(t)> and g2(t) for w2=R
 
for k1 = 1:length(T1)
        
    k1
    
    % <n2(t)>
 
    n2_f1mean_tmp(:,:,k1) = n2_f1(:,:,k1)...
          * rho_ss(:,:,1)/trace(rho_ss(:,:,1));
    n2_f1mean(k1) = trace(n2_f1mean_tmp(:,:,k1));
    S_n2_f1(k1) = gamma_sens/(2*pi*epsilon2^2) * n2_f1mean(k1);
    
    
    % g2(t)

    corr1_tmp(:,:,k1) = s1' * s1 * n2_f1(:,:,k1)...
        * rho_ss(:,:,1)/trace(rho_ss(:,:,1));
    corr1(k1) = trace(corr1_tmp(:,:,k1));
    g2(k1) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2*corr1(k1);
    g2_norm(k1) = g2(k1) / (S_n1(1)*S_n2_f1(k1));
end
  
%% Calculation of <n2(t)> and g2(t) for w2=-R

for k2 = 1:length(T2)
        
    k2
    
    % <n2(t)>

    n2_f2mean_tmp(:,:,k2) = n2_f2(:,:,k2)...
          * rho_ss(:,:,2)/trace(rho_ss(:,:,2));
    n2_f2mean(k2) = trace(n2_f2mean_tmp(:,:,k2));
    S_n2_f2(k2) = gamma_sens/(2*pi*epsilon2^2) * n2_f2mean(k2);
    
    % g2(t)
    
    corr2_tmp(:,:,k2) = s1' * s1 * n2_f2(:,:,k2)...
        * rho_ss(:,:,2)/trace(rho_ss(:,:,2));
    corr2(k2) = trace(corr2_tmp(:,:,k2));
    g2_bis(k2) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2*corr2(k2);
    g2_norm_bis(k2) = g2_bis(k2) / (S_n1(2)*S_n2_f2(k2));
end

%% Plot of <n(t)>

figure; 

colorspec = {[0 0.5 0]; [1 0.8 0]};  % green,yellow

plot(T1,S_n2_f1(:),'Color', colorspec{1});
hold on;
plot(T2,S_n2_f2(:),'Color', colorspec{2});
    
legend('\omega_1=R_2^-','\omega_1=-R',...
        'fontsize', 12,'Location','best');

xlabel('\tau','color','k','fontsize', 20);
ylabel('<n_2(\tau)>','color','k','fontsize', 20);

%% Plot of g2

figure;

u = 1;

line(T1,g2_norm(:),'Color', colorspec{1});
hold on;
line(T2,g2_norm_bis(:),'Color', colorspec{2});

adiff = gca;                    % current axes
adiff_pos = adiff.Position; % position of first axes
axis([0 60 0 6]);
    
legend('\omega_1=R_2^-','\omega_1=R',...
        'fontsize', 12,'Location','best');

xlabel('g \tau','color','k','fontsize', 20);
ylabel('g^{(2)}_{\gamma_2}(\tau)','color','k','fontsize', 20);

hold on;

au = axes('Position',adiff_pos,...
            'XAxisLocation','top','YAxisLocation','right',...
            'YTick',u,'Color','none');
axis([0 60 0 6]);

toc


%% to verify that <n2(t)> = <n2(t)_dag>

for j= 1:length(n2_f1mean)
    
    if n2_f1mean(j) ~= n2_f1mean(j).';    %to write tilde:Alt+5
       disp(['The element ', num2str(j),' is different']) 
    else disp('<n2(t)> = <n2(t)_dag>')
    end
end

%% to check the period of my odd oscillations
% I can calculate the Fourier transform of g2_norm

mm1 = length(g2_norm);          % Window length
nn1 = pow2(nextpow2(mm1));      % Transform length
mm2 = length(g2_norm_bis);          
nn2 = pow2(nextpow2(mm2));

FT1 = fft(g2_norm,nn1);
FT2 = fft(g2_norm_bis,nn2);

DeltaT1 = T1(2) - T1(1);
fs1 = 1/DeltaT1;                  % Sample frequency
x1 = (0:nn1-1)*(fs1/nn1);         % Frequency range

DeltaT2 = T2(2) - T2(1);
fs2 = 1/DeltaT2;
x2 = (0:nn2-1)*(fs2/nn2);

power1 = FT1.*conj(FT1)/nn1;   % Power of the DFT
power2 = FT2.*conj(FT2)/nn2;

figure;

plot(x1,FT1(:),'Color', colorspec{1})
hold on;
plot(x2,FT2(:),'Color', colorspec{2})

legend('configuration ii)','configuration iv)',...
        'fontsize', 12,'Location','best');
xlabel('\omega','color','k','fontsize', 20);
ylabel('g^{(2)}_{\gamma_2}(\omega)','color','k','fontsize', 20);

figure;

plot(x1,power1(:),'Color', 'b')
hold on;
plot(x2,power2(:),'Color', 'r')

legend('configuration ii)','configuration iv)',...
        'fontsize', 12,'Location','best');
xlabel('\omega','color','k','fontsize', 20);
ylabel('|FT|','color','k','fontsize', 20);

%% NO

% newg2_1 = ifft(FT1(5:end));
% newT1 = linspace(0,60, length(FT1(5:end)));
% 
% figure;
% plot(newT1, newg2_1(:));






