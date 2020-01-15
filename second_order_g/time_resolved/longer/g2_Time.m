%% Calculation of the 2 steady states for configurations ii and iv

tic

clear;

rho_ss_tmp = func_ss;

rho_ss = reshape(rho_ss_tmp,[sqrt(size(rho_ss_tmp,1)),...
        sqrt(size(rho_ss_tmp,1)), size(rho_ss_tmp,2)]);

n_configurations = 2;

rho_ss_norm = zeros(size(rho_ss,1),size(rho_ss,2),n_configurations);
rho_ss_norm_vec = reshape(rho_ss_norm,[],n_configurations);

for j=1:n_configurations
    
    rho_ss_norm(:,:,j) = rho_ss(:,:,j)/trace(rho_ss(:,:,j));
end

% % disp('L_sens=0');
% tmp=rho_ss(:,:,1)-rho_ss(:,:,1)'; 
% sum(abs(tmp(:))); disp('sum(abs(tmp))='); disp(sum(abs(tmp(:))));
% max(abs(tmp(:))); disp('max(abs(tmp))='); disp(max(abs(tmp(:))));
% % figure; pcolor(abs(tmp));

load parameters.mat;
    
%% Positive times

%Calculation of <n1(0)> 

load my_operators.mat n ntrunc s1 s2;

    % I'm loading only some of the operators stored in the file
    % my_operators.mat

n1_mean_tmp = zeros(size(rho_ss,1),size(rho_ss,2),n_configurations);
n1_mean = zeros(1,n_configurations);
S_n1 = zeros(1,n_configurations);

for j = 1:n_configurations
    
n1_mean_tmp(:,:,j) = s1' * s1* rho_ss_norm(:,:,j);
n1_mean(j) = trace(n1_mean_tmp(:,:,j));
S_n1(j) = gamma_sens/(2*pi*epsilon1^2) .* n1_mean(j);

end

%% Dynamics of n2(t) with initial condition n2(0)

% ii)

n2 = s2'* s2;
n2_in_vector = reshape(n2,[],1);

L_ii = L_gen(R2_minus,R);
Ldag_ii = Ldag_gen(R2_minus,R);
defun_ii = @(t,Y) Ldag_ii * Y;

[T1,n2_T1] = ode45(defun_ii,[0 60],n2_in_vector);

n2_T1_bis = reshape(n2_T1,[size(n2_T1,1),sqrt(size(n2_T1,2)),...
    sqrt(size(n2_T1,2))]);    % it's a solid figure tx56x56
     
for time = 1:size(n2_T1_bis,1)
   
    n2_ii = permute(n2_T1_bis,[2 3 1]);  % it's a solid figure 56x56xt
end

% iv)

L_iv = L_gen(R,-R);
Ldag_iv = Ldag_gen(R,-R);
defun_iv = @(t,Y) Ldag_iv * Y;

[T2,n2_T2] = ode45(defun_iv,[0 60],n2_in_vector);

n2_T2_bis = reshape(n2_T2,[size(n2_T2,1),sqrt(size(n2_T2,2)),...
    sqrt(size(n2_T2,2))]);
     
for time2 = 1:size(n2_T2_bis,1)
   
    n2_iv = permute(n2_T2_bis,[2 3 1]);
end

%% Calculation of <n2(t)> and g2(t) 
 
% ii)

n2_f1mean_tmp = zeros(size(n2_ii,1),size(n2_ii,2),length(T1));
n2_f1mean = zeros(1,length(T1));
S_n2_ii = zeros(1,length(T1));

corr1_tmp = zeros(size(n2_ii,1),size(n2_ii,2),length(T1));
corr1 = zeros(1,length(T1));
g2_ii = zeros(1,length(T1));
g2_ii_norm = zeros(1,length(T1));

for k1 = 1:length(T1)
        
    k1
    
    % <n2(t)>
 
    n2_f1mean_tmp(:,:,k1) = n2_ii(:,:,k1) * rho_ss_norm(:,:,1);
    n2_f1mean(k1) = trace(n2_f1mean_tmp(:,:,k1));
    S_n2_ii(k1) = gamma_sens/(2*pi*epsilon2^2) .* n2_f1mean(k1);
    
    % g2(t)

    %corr1_tmp(:,:,k1) = (s1' * s1) * n2_ii(:,:,k1) * rho_ss_norm(:,:,1);
    corr1_tmp(:,:,k1) = s1' * n2_ii(:,:,k1) * s1 * rho_ss_norm(:,:,1);
    corr1(k1) = trace(corr1_tmp(:,:,k1));
    g2_ii(k1) = (gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2.*corr1(k1);
    g2_ii_norm(k1) = g2_ii(k1) / (S_n1(1)*S_n2_ii(k1));
end
  
%  iv)

n2_f2mean_tmp = zeros(size(n2_iv,1),size(n2_iv,2),length(T2));
n2_f2mean = zeros(1,length(T2));
S_n2_iv = zeros(1,length(T2));

corr2_tmp = zeros(size(n2_iv,1),size(n2_iv,2),length(T2));
corr2 = zeros(1,length(T2));
g2_iv = zeros(1,length(T2));
g2_iv_norm = zeros(1,length(T2));

for k2 = 1:length(T2)
        
    k2
    
    % <n2(t)>

    n2_f2mean_tmp(:,:,k2) = n2_iv(:,:,k2) * rho_ss_norm(:,:,2);
    n2_f2mean(k2) = trace(n2_f2mean_tmp(:,:,k2));
    S_n2_iv(k2) = gamma_sens/(2*pi*epsilon2^2) .* n2_f2mean(k2);
    
    % g2(t)
    
    %corr2_tmp(:,:,k2) = (s1' * s1) * n2_iv(:,:,k2) * rho_ss_norm(:,:,2);
    corr2_tmp(:,:,k2) = s1' * n2_iv(:,:,k2) * s1 * rho_ss_norm(:,:,2);
    corr2(k2) = trace(corr2_tmp(:,:,k2));
    g2_iv(k2) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2.*corr2(k2);
    g2_iv_norm(k2) = g2_iv(k2) / (S_n1(2)*S_n2_iv(k2));
end

%% NEGATIVE TIMES

%Calculation of <n2(0)>

n2_mean_tmp = zeros(size(rho_ss,1),size(rho_ss,2),n_configurations);
n2_mean = zeros(1,n_configurations);
S_n2 = zeros(1,n_configurations);

for j = 1:n_configurations
    
n2_mean_tmp(:,:,j) = s2' * s2* rho_ss_norm(:,:,j);
n2_mean(j) = trace(n2_mean_tmp(:,:,j));
S_n2(j) = gamma_sens/(2*pi*epsilon1^2) .* n2_mean(j);

end

%% Dynamics of n1(t) with initial condition n1(0)

% ii)

n1 = s1'* s1;
n1_in_vector = reshape(n1,[],1);

[T1_neg,n1_T1] = ode45(defun_ii,[0 60],n1_in_vector);

n1_T1_bis = reshape(n1_T1,[size(n1_T1,1),sqrt(size(n1_T1,2)),...
    sqrt(size(n1_T1,2))]);    % it's a solid figure tx56x56
     
for time_neg = 1:size(n1_T1_bis,1)
   
    n1_ii = permute(n1_T1_bis,[2 3 1]);  % it's a solid figure 56x56xt
end

% iv)

[T2_neg,n1_T2] = ode45(defun_iv,[0 60],n1_in_vector);

n1_T2_bis = reshape(n1_T2,[size(n1_T2,1),sqrt(size(n1_T2,2)),...
    sqrt(size(n1_T2,2))]);
     
for time2_neg = 1:size(n1_T2_bis,1)
   
    n1_iv = permute(n1_T2_bis,[2 3 1]);
end

%% Calculation of <n1(t)> and g2(t) 
 
% ii)

n1_f1mean_tmp = zeros(size(n1_ii,1),size(n1_ii,2),length(T1_neg));
n1_f1mean = zeros(1,length(T1_neg));
S_n1_ii = zeros(1,length(T1_neg));

corr1n_tmp = zeros(size(n2_ii,1),size(n2_ii,2),length(T1_neg));
corr1n = zeros(1,length(T1_neg));
g2n_ii = zeros(1,length(T1_neg));
g2n_ii_norm = zeros(1,length(T1_neg));

for k1 = 1:length(T1_neg)
        
    k1
    
    % <n2(t)>
 
    n1_f1mean_tmp(:,:,k1) = n1_ii(:,:,k1) * rho_ss_norm(:,:,1);
    n1_f1mean(k1) = trace(n1_f1mean_tmp(:,:,k1));
    S_n1_ii(k1) = gamma_sens/(2*pi*epsilon2^2) .* n1_f1mean(k1);
    
    % g2(t)

    %corr1n_tmp(:,:,k1) = n1_ii(:,:,k1) * (s2' * s2) * rho_ss_norm(:,:,1);
    corr1n_tmp(:,:,k1) = s2' * n1_ii(:,:,k1) * s2 * rho_ss_norm(:,:,1);
    corr1n(k1) = trace(corr1n_tmp(:,:,k1));
    g2n_ii(k1) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2.*corr1n(k1);
    g2n_ii_norm(k1) = g2n_ii(k1) / (S_n2(1).*S_n1_ii(k1));
end
  
% iv)

n1_f2mean_tmp = zeros(size(n2_iv,1),size(n2_iv,2),length(T2_neg));
n1_f2mean = zeros(1,length(T2_neg));
S_n1_iv = zeros(1,length(T2_neg));

corr2n_tmp = zeros(size(n2_iv,1),size(n2_iv,2),length(T2_neg));
corr2n = zeros(1,length(T2_neg));
g2n_iv = zeros(1,length(T2_neg));
g2n_iv_norm = zeros(1,length(T2_neg));

for k2 = 1:length(T2_neg)
        
    k2
    
    % <n2(t)>

    n1_f2mean_tmp(:,:,k2) = n1_iv(:,:,k2) * rho_ss_norm(:,:,2);
    n1_f2mean(k2) = trace(n1_f2mean_tmp(:,:,k2));
    S_n1_iv(k2) = gamma_sens/(2*pi*epsilon2^2) .* n1_f2mean(k2);
    
    % g2(t)
  
    %corr2n_tmp(:,:,k2) = n1_iv(:,:,k2) * (s2' * s2) * rho_ss_norm(:,:,2);
    corr2n_tmp(:,:,k2) = s2' * n1_iv(:,:,k2) * s2 * rho_ss_norm(:,:,2);
    corr2n(k2) = trace(corr2n_tmp(:,:,k2));
    g2n_iv(k2) =(gamma_sens^2)/(2*pi*epsilon1*epsilon2)^2.*corr2n(k2);
    g2n_iv_norm(k2) = g2n_iv(k2) / (S_n2(2).*S_n1_iv(k2));
end

%% check Tr{A*L*rho} = Tr{Ldag*A*rho}

% positive times

tr_n2L_ii = trace(rho_ss_norm_vec(:,1) * n2_in_vector' * L_ii);
tr_n2Ldag_ii = trace(Ldag_ii * n2_in_vector * rho_ss_norm_vec(:,1)');

tr_n2L_iv = trace(rho_ss_norm_vec(:,1) * n2_in_vector' * L_iv);
tr_n2Ldag_iv = trace(Ldag_iv * n2_in_vector * rho_ss_norm_vec(:,1)');

% negative times

tr_n1L_ii = trace(rho_ss_norm_vec(:,1) * n1_in_vector' * L_ii);
tr_n1Ldag_ii = trace(Ldag_ii * n1_in_vector * rho_ss_norm_vec(:,1)');

tr_n1L_iv = trace(rho_ss_norm_vec(:,1) * n1_in_vector' * L_iv);
tr_n1Ldag_iv = trace(Ldag_iv * n1_in_vector * rho_ss_norm_vec(:,1)');

% configuration ii

if tr_n2L_ii==tr_n2Ldag_ii & tr_n1L_ii==tr_n1Ldag_ii 
    disp('Tr{n*L_ii*rho} = Tr{Ldag_ii*n*rho}')
else disp('Tr{n*L_ii*rho} is different from Tr{Ldag_ii*n*rho}')
end

% configuration iv

if tr_n2L_iv==tr_n2Ldag_iv & tr_n1L_iv==tr_n1Ldag_iv 
    disp('Tr{n*L_iv*rho} = Tr{Ldag_iv*n*rho}')
else disp('Tr{n*L_iv*rho} is different from Tr{Ldag_iv*n*rho}')
end

%% Plot of <n(t)>
% positive and negative times 

figure; 

colorspec = {[0 0.5 0]; [1 0.8 0]};  % green,yellow

Tneg1 = linspace(-60,0,length(T1_neg));
Tneg2 = linspace(-60,0,length(T2_neg));

plot(T1,S_n2_ii(:),'o','Color', colorspec{1});
hold on;
plot(T2,S_n2_iv(:),'Color', colorspec{2});
hold on;

legend('ii)','iv)','fontsize', 12,'Location','best');

plot(Tneg1,S_n1_ii(:),'Color', colorspec{1});
hold on;
plot(Tneg2,S_n1_iv(:),'Color', colorspec{2});

xlabel('t','color','k','fontsize', 20);
ylabel('<n(t)>','color','k','fontsize', 20);

%% Plot of g2
% positive and negative times

figure;

u = 1;

line(T1,g2_ii_norm(:),'Color', colorspec{1});
hold on;
line(T2,g2_iv_norm(:),'Color', colorspec{2});
hold on;

legend('ii)','iv)','fontsize', 12,'Location','northwest');
    
%line(Tneg1,flip(g2n_ii_norm(:)),'Color', colorspec{1});
line(Tneg1,g2n_ii_norm(end:-1:1),'Color', colorspec{1});
hold on;
%line(Tneg2,flip(g2n_iv_norm(:)),'Color', colorspec{2});
line(Tneg2,g2n_iv_norm(end:-1:1),'Color', colorspec{2});

toc

adiff = gca;                    % current axes
adiff_pos = adiff.Position;     % position of first axes
axis([-60 60 -0.1 3.4]);

xlabel('g \tau','color','k','fontsize', 20);
ylabel('g^{(2)}_{\gamma_2}(\tau)','color','k','fontsize', 20);

hold on;

au = axes('Position',adiff_pos,...
            'XAxisLocation','top','YAxisLocation','right',...
            'YTick',u,'Color','none');
axis([-60 60 -0.1 3.4]);

%toc