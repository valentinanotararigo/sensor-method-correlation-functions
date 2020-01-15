%% Calculation of the 2 steady states for configurations ii and iv

tic

clear;

rho_ss_tmp = func_ss;

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)), size(rho_ss_tmp,2)]);
    
n_configurations = 2;

for j=1:n_configurations
    
    rho_ss(:,:,j) = rho_ss(:,:,j)/trace(rho_ss(:,:,j));
end

load parameters.mat;
    
%% Calculation of <n2(0)>

load my_operators.mat n ntrunc s1 s2;

    % I'm loading only some of the operators stored in the file
    % my_operators.mat

for j = 1:n_configurations
    
n2_mean_tmp(:,:,j) = s2' * s2 * rho_ss(:,:,j)/trace(rho_ss(:,:,j));
n2_mean(j) = trace(n2_mean_tmp(:,:,j));
S_n2(j) = gamma_sens/(2*pi*epsilon1^2) * n2_mean(j);

end

%% Dynamics of n1(t) with initial condition n1(0)

% w1=R2^-

n1 = s1'* s1;
n1_in_vector = reshape(n1,[],1);

[T1,n1_T1] = ode45(@func_n2_T1,[0 60],n1_in_vector);

n1_T1_bis = reshape(n1_T1,[size(n1_T1,1),sqrt(size(n1_T1,2)),...
    sqrt(size(n1_T1,2))]);    % it's a solid figure tx56x56
     
for time = 1:size(n1_T1_bis,1)
   
    n1_f1 = permute(n1_T1_bis,[2 3 1]);  % it's a solid figure 56x56xt
end

% w1=R

[T2,n1_T2] = ode45(@func_n2_T2,[0 60],n1_in_vector);

n1_T2_bis = reshape(n1_T2,[size(n1_T2,1),sqrt(size(n1_T2,2)),...
    sqrt(size(n1_T2,2))]);
     
for time2 = 1:size(n1_T2_bis,1)
   
    n1_f2 = permute(n1_T2_bis,[2 3 1]);
end

%% Calculation of <n1(t)> and g2(t) for configuration ii
 
for k1 = 1:length(T1)
    
    % <n2(t)>
 
    n1_f1mean_tmp(:,:,k1) = n1_f1(:,:,k1) * rho_ss(:,:,1);
    n1_f1mean(k1) = trace(n1_f1mean_tmp(:,:,k1));
    S_n1_f1(k1) = gamma_sens/(2*pi*epsilon2^2) * n1_f1mean(k1);
    
    % disp('L_sens=0');
    tmp_1(:,:,k1)=n1_f1(:,:,k1)-n1_f1(:,:,k1)'; 
    tmp_n1(:,k1) = reshape(tmp_1(:,:,k1), size(tmp_1,1)^2,[]);
    sum(abs(tmp_n1(:,k1))); 
    max(abs(tmp_n1(:,k1))); 
    if sum(abs(tmp_n1(:,k1))) > 10^-13
        disp(['At T1=',num2str(k1),... 
            ' sum(abs(tmp)) = ', num2str(sum(abs(tmp_n1(:,k1))))]);
        disp(['And max(abs(tmp)) = ', num2str(max(abs(tmp_n1(:,k1))))]);
    end
    % figure; pcolor(abs(tmp));
    
    
    % g2(t)

    corr1n_tmp(:,:,k1) = n1_f1(:,:,k1) * s2' * s2 * rho_ss(:,:,1);
    corr1n(k1) = trace(corr1n_tmp(:,:,k1));
    g2n(k1) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2*corr1n(k1);
    g2n_norm(k1) = g2n(k1) / (S_n2(1)*S_n1_f1(k1));
end

%% Calculation of <n1(t)> and g2(t) for configuration iv

for k2 = 1:length(T2)
        
    k2
    
    % <n2(t)>

    n1_f2mean_tmp(:,:,k2) = n1_f2(:,:,k2) * rho_ss(:,:,2);
    n1_f2mean(k2) = trace(n1_f2mean_tmp(:,:,k2));
    S_n1_f2(k2) = gamma_sens/(2*pi*epsilon2^2) * n1_f2mean(k2);
    
    % g2(t)
  
    corr2n_tmp(:,:,k2) = n1_f2(:,:,k2) * s2' * s2 * rho_ss(:,:,2);
    corr2n(k2) = trace(corr2n_tmp(:,:,k2));
    g2n_bis(k2) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2*corr2n(k2);
    g2n_norm_bis(k2) = g2n_bis(k2) / (S_n2(2)*S_n1_f2(k2));
end

%% Plot of <n(t)>

figure; 

colorspec = {[0 0.5 0]; [1 0.8 0]};  % green,yellow

Tneg1 = linspace(-60,0,length(T1));
Tneg2 = linspace(-60,0,length(T2));

plot(Tneg1,S_n1_f1(:),'Color', colorspec{1});
hold on;
plot(Tneg2,S_n1_f2(:),'Color', colorspec{2});
    
legend('\omega_1=R_2^-','\omega_1=-R',...
        'fontsize', 12,'Location','best');

xlabel('\tau','color','k','fontsize', 20);
ylabel('<n_1(\tau)>','color','k','fontsize', 20);

%% Plot of g2

figure;

u = 1;

line(Tneg1,flip(g2n_norm(:)),'Color', colorspec{1});
hold on;
line(Tneg2,flip(g2n_norm_bis(:)),'Color', colorspec{2});

adiff = gca;                    % current axes
adiff_pos = adiff.Position; % position of first axes
axis([-60 0 -1.5 2]);
    
legend('\omega_1=R_2^-','\omega_1=-R',...
        'fontsize', 12,'Location','best');

xlabel('gt','color','k','fontsize', 20);
ylabel('g^{(2)}_{\gamma_2}(t)','color','k','fontsize', 20);

hold on;

au = axes('Position',adiff_pos,...
            'XAxisLocation','top','YAxisLocation','right',...
            'YTick',u,'Color','none');
axis([-60 0 -1.5 2]);

toc