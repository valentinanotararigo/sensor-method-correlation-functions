%% 

tic

Lfromfunc = func_sensor_jk;

% rhotot_in_vector = zeros(size(Lfromfunc,1),1);
% rhotot_in_vector(1) = 1;

%%
g=1;
nloop = 200;
gamma_mode = [0.01; 0.1; 0.5]*g;

rho_ss_tmp = zeros(size(Lfromfunc,1),nloop,length(gamma_mode));

for j=1:length(gamma_mode)
    
    for k = 1:nloop
    
        k
    
        if k ==1
            [A,D] = eigs(Lfromfunc(:,:,j,k),1,0);    
            %eigs(A,k,sigma) returns k eigenvalues with value=sigma
        else
            [A,D] = eigs(Lfromfunc(:,:,j,k),1,0,opts);
            %eigs(A,K,sigma,opts) specifies an options structure. 
        end
    
        rho_ss_tmp(:,j,k) = A;
        opts.v0 = A;           %opts.v0=starting vector
    end
end

%%
n=6;
ntrunc = [1,n-1,1];

s = func_operators(ntrunc,3);

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)),size(rho_ss_tmp,2),size(rho_ss_tmp,2)]);

w_mode = 0;
g = 1;

gamma_atom = 0.01*g;
gamma_sens = 0.001*g;
epsilon = 0.00001;

diff = linspace(-3.5,3.5,nloop);   %diff=(w_sens-w_mode)/g
w_sens = g.*diff + w_mode;
%%
rho_ss_norm = zeros(length(gamma_mode),nloop);
rho_corr_tmp = zeros(size(rho_ss,1),...
    size(rho_ss,2),length(gamma_mode),nloop);
rho_corr = zeros(length(gamma_mode),nloop);
S = zeros(length(gamma_mode),nloop);

for j = 1:length(gamma_mode)
    
    for k = 1:nloop
    
        rho_ss_norm(j,k) = trace(rho_ss(:,:,j,k));
        rho_corr_tmp(:,:,j,k) = s' * s * rho_ss(:,:,j,k)/rho_ss_norm(j,k);
        rho_corr(j,k) = trace(rho_corr_tmp(:,:,j,k));
        S(j,k) = gamma_sens/(2*pi*epsilon^2) * rho_corr(j,k);
    
    end
end
%%

u = [-3.16 -2.415 -1 -0.42 -0.31 0.3 0.41 1 2.42 3.15];

figure;
hold on;

colorspec = {[0 0 1]; [1 0 0]; [1 0.8 0]};  % blue,red,yellow

for j = 1:length(gamma_mode)

    line(diff,log(S(j,:)),'Color', colorspec{j},'LineWidth',1);
    
    adiff = gca; % current axes
    adiff_pos = adiff.Position; % position of first axes
    axis([-3.5 3.5 -11.5 1.5]);
    set(gca,'fontsize', 15);
    legend('\gamma_a=0.01g','\gamma_a=0.1g',...
        '\gamma_a=0.5g','Location','best');
end

% title('Power spectrum',...
%     'color','k','fontsize', 18,'fontname','helvetica',...
%     'fontunits','normalized','fontweight','normal');

xlabel('(\omega_1-\omega_a)/g','color','k','fontsize', 20);
ylabel('S^{(1)}_{\Gamma_0}(\omega_1)','color','k','fontsize', 20);

hold on;

%plot(u,ones(1,length(u)));
au = axes('Position',adiff_pos,'XAxisLocation','top',...
    'YAxisLocation','right','XTick',u,'Ytick',[],...
    'XTickLabel',{'-R_3^+','-R_2^+','-R','-R_2^-','-R_3^-',...
    'R_3^-','R_2^-','R','R_2^+','R_3^+'},...
    'XGrid','on','Color','none');
set(au,'fontsize', 15);
line(u,ones(1,length(u)),'Color','none','Parent',au);
axis([-3.5 3.5 -11.5 1.5]);

toc