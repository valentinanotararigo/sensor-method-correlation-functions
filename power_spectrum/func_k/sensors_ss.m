%% 

tic

Lfromfunc = func_sensors_k;

% rhotot_in_vector = zeros(size(Lfromfunc,1),1);
% rhotot_in_vector(1) = 1;

%%

nloop = 200;

A = zeros(size(Lfromfunc,1),size(Lfromfunc,1),nloop);
D = zeros(size(Lfromfunc,1),size(Lfromfunc,1),nloop);

x = zeros(size(Lfromfunc,1),nloop);
x_ascendent = zeros(size(Lfromfunc,1),nloop);
clear lg rho_ss_tmp

for k = 1:nloop
    
    k
    
    [A(:,:,k),D(:,:,k)] = eig(Lfromfunc(:,:,k));
    
    x(:,k) = diag(D(:,:,k));
    
        % x=diag(A) returns a column vector of the main diagonal elements 
        % of A. In my case x is a matrix, where each column contains the
        % eigenvalues of L, corresponding to a different k.
        
    x_ascendent(:,k) = sort(x(:,k));
    
        % B = sort(A) sorts the elements of A in ascending order along the 
        % first array dimension whose size does not equal 1. In my case A 
        % is a matrix, so sort(A) put in ascending order the elements of
        % each column.
        
    lg(:,k) = abs(x(:,k)) < eps(10^10);
    
        % eps returns the distance from 1.0 to the next largest 
        % double-precision number, that is, eps = 2^-52.
        % eps(x) returns the positive distance from abs(x) to the next 
        % largest floating-point number of the same precision as x
        
    rho_ss_tmp = A(:,lg);
end

%%
n=6;
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

rho_corr_tmp = zeros(sqrt(size(rho_ss_tmp,1)),...
    sqrt(size(rho_ss_tmp,1)),size(rho_ss_tmp,2));
rho_corr = zeros(nloop,1);
S = zeros(nloop,1);

for k = 1:nloop
    
    rho_corr_tmp(:,:,k) = s' * s * rho_ss(:,:,k);
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