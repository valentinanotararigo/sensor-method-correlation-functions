%% 

tic

Lfromfunc = func_sensors_k;

% rhotot_in_vector = zeros(size(Lfromfunc,1),1);
% rhotot_in_vector(1) = 1;

%%

nloop = 100;

%A = zeros(size(Lfromfunc,1),size(Lfromfunc,1),nloop);
%D = zeros(size(Lfromfunc,1),size(Lfromfunc,1),nloop);
rho_ss_tmp = zeros(size(Lfromfunc,1),nloop);

for k = 1:nloop
    
    k
    
    if k ==1
        [A,D] = eigs(Lfromfunc(:,:,k),1,0);    
        %eigs(A,k,sigma) returns k eigenvalues with value=sigma
    else
        [A,D] = eigs(Lfromfunc(:,:,k),1,0,opts);
        %eigs(A,K,sigma,opts) specifies an options structure. 
    end
    
    rho_ss_tmp(:,k) = A;
    opts.v0 = A;           %opts.v0=starting vector
end

%%
n=3;
ntrunc = [1,n-1,1];

s = func_operators(ntrunc,3);

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)),size(rho_ss_tmp,2)]);

w_mode = 0;
g = 1;

gamma_atom = 0.01*g;
gamma_mode = 0.5*g;
gamma_sens = 0.001*g;
epsilon = 0.00001;

diff = linspace(-3,3,nloop);   %diff=(w_sens-w_mode)/g
w_sens = g.*diff + w_mode;

rho_ss_norm = zeros(1,nloop);
rho_corr_tmp = zeros(sqrt(size(rho_ss_tmp,1)),...
    sqrt(size(rho_ss_tmp,1)),size(rho_ss_tmp,2));
rho_corr = zeros(nloop,1);
S = zeros(nloop,1);

for k = 1:nloop
    
    rho_ss_norm(k) = trace(rho_ss(:,:,k));
    rho_corr_tmp(:,:,k) = s' * s * rho_ss(:,:,k)/rho_ss_norm(k);
    rho_corr(k,1) = trace(rho_corr_tmp(:,:,k));
    S(k,1) = gamma_sens/(2*pi*epsilon^2) * rho_corr(k,1);
    
end

%%

figure;
semilogy(diff,S);

legend('\gamma_{mode}=0.5g','Location','bestoutside');

title('Power spectrum',...
    'color','k','fontsize', 18,'fontname','helvetica',...
    'fontunits','normalized','fontweight','normal');

xlabel('(\omega_{sens}-\omega_{mode})/g','color','k','fontsize', 12);
ylabel('S','color','k','fontsize', 12);

toc