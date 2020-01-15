% Function to get the Liouvillian and to solve the Linblad master equation 
% of a two level system coupled to a resonant mode, considering also the 
% interaction with sensors.

function dY = func_sensors_k(T,Y)

persistent L

if nargin==0

n=3;   % to calculate sigma_z I need a_tmp
ntrunc = [1,n-1,1];

sigma_m = func_operators(ntrunc,1);
a = func_operators(ntrunc,2);
s = func_operators(ntrunc,3);

sigma_p = sigma_m';
adag = a';     % a is the transpose of adag

w_mode = 0;                
ht=1;
g = 1;

H_JC = ht * (sigma_p*a + sigma_m*adag);

gamma_atom = 0.01*g;
gamma_mode = 0.5*g;
gamma_pump = gamma_atom;
gamma_sens = 0.001*g;
epsilon = 0.00001;

H_int = ht * epsilon * (a*s'+a'*s);

H_inv = H_JC + H_int;

L_atom = gamma_atom * (kron(sigma_m,sigma_m) ...
         - 0.5 * kron(sigma_p*sigma_m,eye(length(H_inv))) ...
         - 0.5 * kron(eye(length(H_inv)),sigma_p*sigma_m));
     
L_mode = gamma_mode * (kron(a,a) ...
        - 0.5 * kron(adag*a,eye(length(H_inv))) ...
        - 0.5 * kron(eye(length(H_inv)),adag*a));
    
L_pump = gamma_pump * (kron(sigma_p,sigma_p) ...
         - 0.5 * kron(sigma_m*sigma_p,eye(length(H_inv)))...
         - 0.5 * kron(eye(length(H_inv)),sigma_m*sigma_p));
     
L_inv = L_atom + L_mode + L_pump;

nloop = 100;

diff = linspace(-3,3,nloop);   %diff=(w_sens-w_mode)/g
w_sens = g.*diff + w_mode;

for k = 1:nloop
    
    H_sens = ht * w_sens(k) * s' * s;
    
    H = H_inv + H_sens;
    
    LH = -1i * (kron(H.',eye(length(H)))-kron(eye(length(H)),H));
    L_sens = gamma_sens * (kron(s,s) ...
         - 0.5 * kron(s'*s,eye(length(H_sens))) ...
         - 0.5 * kron(eye(length(H_sens)),s'*s));

    L(:,:,k) = LH + L_inv + L_sens;

end

dY = L; return

end

dY = L*Y;
