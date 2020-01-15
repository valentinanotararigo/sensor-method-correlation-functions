function dY = func_sensors(T,Y)

n=3;   % to calculate sigma_z I need a_tmp
ntrunc = [1,n-1,1];

sigma_m = func_operators(ntrunc,1);
a = func_operators(ntrunc,2);
s = func_operators(ntrunc,3);

sigma_p = sigma_m';
adag = a';     % a is the transpose of adag

l = 1:n;
a_tmp = sparse(1:n,2:n+1,sqrt(l));
a_tmp(:,size(a_tmp,2))=[]; %to eliminate the last column
sigma_z_tmp = [1 0; 0 -1];
sigma_z = kron( kron(sigma_z_tmp,eye(length(a_tmp))), ...
    eye(length(sigma_z_tmp)));

%w_mode = 1000;            % cavity frequency given in cm^-1 (infrared)
w_mode = 2;
w_atom = w_mode;                    
ht=1;
%g = 267.1;  % like in my notes: with this value RWA is valid
g = 0.5;

H_atom = 0.5 * ht * w_atom * sigma_z; 
H_mode = ht * w_mode * adag * a;
H_JC = ht * g * (sigma_p*a + sigma_m*adag);
%H_JC = ht * (sigma_p*a + sigma_m*adag);

gamma_atom = 0.01 * g;
gamma_mode = 0.5 * g;
gamma_pump = 0;%gamma_atom;
gamma_sens = 0.1 * g;
epsilon = (sqrt(gamma_sens*gamma_mode)/2) /10;

H_int = ht * epsilon * (a*s'+a'*s);

H_inv = H_atom + H_mode + H_JC + H_int;

L_atom = gamma_atom * (kron(sigma_m,sigma_m) ...
         - 0.5 * kron(sigma_p*sigma_m,eye(length(H_atom))) ...
         - 0.5 * kron(eye(length(H_atom)),sigma_p*sigma_m));
     
L_mode = gamma_mode * (kron(a,a) ...
        - 0.5 * kron(adag*a,eye(length(H_mode))) ...
        - 0.5 * kron(eye(length(H_mode)),adag*a));
    
L_pump = gamma_pump * (kron(sigma_p,sigma_p) ...
         - 0.5 * kron(sigma_m*sigma_p,eye(length(H_atom)))...
         - 0.5 * kron(eye(length(H_atom)),sigma_m*sigma_p));
     
L_inv = L_atom + L_mode + L_pump;

%diff = linspace(-3,3,nloop);   %diff=(w_sens-w_mode)/g
diff = 1;
w_sens = g*diff + w_mode;

H_sens = ht * w_sens * s' * s;
%H_sens = ht * diff(k) * s' * s;
    
H = H_inv + H_sens;
    
LH = -1i * (kron(H.',eye(length(H)))-kron(eye(length(H)),H));
    
L_sens = gamma_sens * (kron(s,s) ...
         - 0.5 * kron(s'*s,eye(length(H_sens))) ...
         - 0.5 * kron(eye(length(H_sens)),s'*s));

L = LH + L_inv + L_sens;

dY = L*Y;
