
% steady state from

function rss_norm = rss_dyn(L,tspan,r_in) 

defunc = @(x,z) L * z;              
[t,y] = ode45(defunc,tspan,r_in);  
rho_ss = y(length(t),:).';     % vec
rss_matr = reshape(rho_ss,sqrt(length(rho_ss)),sqrt(length(rho_ss)));
rss_norm = rss_matr / trace(rss_matr);   % matrix

