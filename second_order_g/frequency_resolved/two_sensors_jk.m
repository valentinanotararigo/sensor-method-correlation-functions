%% 

tic

parameters;
my_operators;

rho_ss_tmp = func_2sensors_jk;

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)),size(rho_ss_tmp,2),size(rho_ss_tmp,3)]);
    
%%

rho_ss_norm = zeros(size(rho_ss,1),size(rho_ss,2),length(gamma_mode),nloop);
rho_corr1_tmp = zeros(size(rho_ss,1),size(rho_ss,2),...
                length(gamma_mode),nloop);
rho_corr1 = zeros(length(gamma_mode),nloop);
n1 = zeros(length(gamma_mode),nloop);
rho_corr2_tmp = zeros(size(rho_ss,1),size(rho_ss,2),...
                length(gamma_mode),nloop);
rho_corr2 = zeros(length(gamma_mode),nloop);
n2 = zeros(length(gamma_mode),nloop);
rho_2ndorder_tmp = zeros(size(rho_ss,1),size(rho_ss,2),...
                   length(gamma_mode),nloop);
rho_2ndorder = zeros(length(gamma_mode),nloop);
g2 = zeros(length(gamma_mode),nloop);
g2_norm = zeros(length(gamma_mode),nloop);

for j = 1:length(gamma_mode)
    
    for k = 1:nloop
        
        rho_ss_norm(:,:,j,k) = rho_ss(:,:,j,k) / trace(rho_ss(:,:,j,k));
        
        rho_corr1_tmp(:,:,j,k) = s1' * s1 * rho_ss_norm(:,:,j,k);
        rho_corr1(j,k) = trace(rho_corr1_tmp(:,:,j,k));
        n1(j,k) = gamma_sens1(j)/(2*pi*epsilon1^2) * rho_corr1(j,k);
        
        rho_corr2_tmp(:,:,j,k) = s2' * s2 * rho_ss_norm(:,:,j,k);
        rho_corr2(j,k) = trace(rho_corr2_tmp(:,:,j,k));
        n2(j,k) = gamma_sens2(j)/(2*pi*epsilon2^2) * rho_corr2(j,k);
    
        rho_2ndorder_tmp(:,:,j,k) = s1' * s1 * s2' * s2...
                                    * rho_ss_norm(:,:,j,k);
        rho_2ndorder(j,k) = trace(rho_2ndorder_tmp(:,:,j,k));
        g2(j,k) = (gamma_sens1(j)*gamma_sens2(j))/...
            (2*pi*epsilon1*epsilon2)^2 * rho_2ndorder(j,k);
        g2_norm(j,k) = g2(j,k) / (n1(j,k)*n2(j,k));
    end
end

%%

u = [-3.16 -2.415 -1 -0.42 -0.31 0.3 0.41 1 2.42 3.15];

figure;
hold on;

colorspec = {[0 0 1]; [1 0 0]; [1 0.8 0]};  % blue,red,yellow

for j = 1:length(gamma_mode)

    line(diff,g2_norm(j,:),'Color', colorspec{j},'LineWidth',1);
    
    adiff = gca; % current axes
    adiff_pos = adiff.Position; % position of first axes
    axis([-3.5 3.5 0 2]);
    set(gca,'fontsize', 15);
    legend('\gamma_a=0.01g','\gamma_a=0.1g',...
        '\gamma_a=0.5g','fontsize', 12,'Location','northeast');
end

xlabel('(\omega_1-\omega_a)/g','color','k','fontsize', 20);
ylabel('g^{(2)}(\omega_1,\omega_2)','color','k','fontsize', 20);

hold on;

au = axes('Position',adiff_pos,'XAxisLocation','top',...
    'YAxisLocation','right','XTick',u,'Ytick',[],...
    'XTickLabel',{'-R_3^+', '-R_2^+', '-R','-R_2^-','-R_3^-',...
    'R_3^-','R_2^-','R','R_2^+','R_3^+'},'XGrid','on','Color','none');
set(au,'fontsize', 15);
line(u,ones(1,length(u)),'Color','none','Parent',au);
axis([-3.5 3.5 0 2]);

toc

