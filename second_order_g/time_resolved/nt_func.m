
% Dynamics of n(t) with initial condition n(0)

function n_time = nt_func(s,Ldag,tau) 

n = s'* s;
n_vec = reshape(n,[],1);
tspan = 0:0.1:tau;

defunc = @(t,Y) Ldag * Y;
[T,nT] = ode45(defunc,tspan,n_vec);

nT_bis = reshape(nT,[size(nT,1),sqrt(size(nT,2)),sqrt(size(nT,2))]);    
% it's a solid figure tx56x56
     
for time = 1:size(nT_bis,1)
    n_time = permute(nT_bis,[2 3 1]);  % it's a solid figure 56x56xt
end

T = tspan;