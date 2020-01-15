%% Calculation of the 2 steady states for configurations ii and iv

tic

clear;

parameters;
my_operators;

L_ii = L_gen(R2_minus,R);
L_iv = L_gen(R,-R);
Ldag_ii = Ldag_gen(R2_minus,R);
Ldag_iv = Ldag_gen(R,-R);

% total steady state

disp('rss');

rss_ii = rss_L(L_ii);                      % matrix
rss_iv = rss_L(L_iv);
% rss_ii_dyn = rss_dyn(L_ii,0:50,rho0);    % matrix
% rss_iv_dyn = rss_dyn(L_iv,0:50,rho0);
% compare_ss(rss_ii,rss_ii_dyn);
% compare_ss(rss_iv,rss_iv_dyn);

% Dynamics of n2(t) and n1(t)

tau = 60;

disp('n2(t)');
n2_ii = nt_func(s2,Ldag_ii,tau);       % tau>0
n2_iv = nt_func(s2,Ldag_iv,tau);
disp('n1(t)');
n1_ii = nt_func(s1,Ldag_ii,tau);       % tau>0
n1_iv = nt_func(s1,Ldag_iv,tau);

%%
% Calculation of g2(t) 

disp('G');
[t,g2p_ii] = g2func(rss_ii,s1,n2_ii,tau);    % tau>0
[t,g2p_iv] = g2func(rss_iv,s1,n2_iv,tau);
[t,g2n_ii] = g2func(rss_ii,s2,n1_ii,tau);    % tau<0
[t,g2n_iv] = g2func(rss_iv,s2,n1_iv,tau);

%% Plot of g2 positive and negative times

figure;

colorspec = {[0 0.5 0]; [1 0.8 0]};  % green,yellow
t_neg = linspace(-60,0,length(t));

plot(t,g2p_ii,'Color', colorspec{1},'LineWidth',1);
hold on;
plot(t,g2p_iv,'Color', colorspec{2},'LineWidth',1);
hold on;

set(gca,'fontsize', 15);
legend('ii)','iv)','fontsize', 12,'Location','northwest');
    
plot(t_neg,g2n_ii(end:-1:1),'Color', colorspec{1},'LineWidth',1);
hold on;
plot(t_neg,g2n_iv(end:-1:1),'Color', colorspec{2},'LineWidth',1);

adiff = gca;                    % current axes
adiff_pos = adiff.Position;     % position of first axes
%axis([-60 60 -0.1 3.5]);     %normal order
axis([-60 60 -2 6]);      % no normal order

xlabel('g \tau','color','k','fontsize', 20);
ylabel('g^{(2)}_{\gamma_2}(\tau)','color','k','fontsize', 20);

hold on;

au = axes('Position',adiff_pos,'XAxisLocation','top',...
            'YAxisLocation','right','XTick',0,'Ytick',1,...
            'XTickLabel',[],'YTickLabel',[],...
            'XGrid','on','YGrid','on','Color','none');
set(au,'fontsize', 15);
%axis([-60 60 -0.1 3.5]);
axis([-60 60 -2 6]);

toc