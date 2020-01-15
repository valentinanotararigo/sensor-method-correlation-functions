
% Calculation of <n(t)> and g(t) 
 
function [tspan,g2_norm] = g2func(yss,s,n_time,tau)

parameters;

n0mean_tmp = s'* s* yss;
n0mean = trace(n0mean_tmp);
S0 = gamma_sens/(2*pi*epsilon^2) * n0mean;

tspan = 0:0.1:tau;

nt_mean_tmp = zeros(size(n_time,1),size(n_time,2),length(tspan));
nt_mean = zeros(length(tspan));
S_t = zeros(length(tspan));
corr_tmp = zeros(size(n_time,1),size(n_time,2),length(tspan));
corr = zeros(length(tspan));
g2 = zeros(length(tspan));
%g2_norm = zeros(length(tspan));

for k = 1:length(tspan)
    
    % <n(t)>
 
    nt_mean_tmp(:,:,k) = n_time(:,:,k) * yss;
    nt_mean(k) = trace(nt_mean_tmp(:,:,k));
    S_t(k) = gamma_sens/(2*pi*epsilon^2) * nt_mean(k);
    
    % g(t)

    corr_tmp(:,:,k) = s' * s * n_time(:,:,k) * yss;
    %corr_tmp(:,:,k) = s' * n_time(:,:,k) * s * yss;
    corr(k) = trace(corr_tmp(:,:,k));
    g2(k) = (gamma_sens*gamma_sens)/((2*pi)^2*epsilon^4) *corr(k);
    g2_norm(k) = g2(k) / (S0*S_t(k));
end