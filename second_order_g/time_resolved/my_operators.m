
% operators

n=5;   
ntrunc = [1,n-1,1,1];  %the two level atom, the mode, 2 sensors

sigma_m = func_operators(ntrunc,1);
sigma_p = sigma_m';
a = func_operators(ntrunc,2);
s1 = func_operators(ntrunc,3);
s2 = func_operators(ntrunc,4);

ground = [1;0];
exc = [0;1];
vacuum = [1;0;0;0;0];

rho_atom = ground * ground';
rho_mode = vacuum * vacuum';
rho_sens = kron(ground,ground) * kron(ground,ground)';  
rho_in = kron(kron(rho_atom,rho_mode),rho_sens);
rho0 = reshape(rho_in,1,[]);

%save my_operators.mat