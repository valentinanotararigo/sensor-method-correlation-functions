% Calculation of the SS with w1=R2^-

% Function to get the Liouvillian and to solve the Linblad master equation 
% of a two level system coupled to a resonant mode, considering also the 
% interaction with 2 sensors.

function dY = func_ss(T,Y)

persistent rho_ss_tmp_func

if nargin==0

load my_operators.mat;

load parameters.mat;

w1 = [R2_minus,R];
w2 = [R,-R];

H_JC = ht * (sigma_p*a + sigma_m*a');
H_int = ht * (epsilon1 * (a*s1'+a'*s1) + epsilon2 * (a*s2'+a'*s2));

for j=1:length(w1)

    H_sens = ht * (w1(j) * s1' * s1 + w2(j) * s2' * s2);
    H = H_JC + H_int + H_sens;
    
    LH = commutator(H);

    L_atom = dissipator(gamma_atom, sigma_m, H);    

    L_pump = dissipator(gamma_pump, sigma_p, H); 

    L_mode = dissipator(gamma_mode, a, H); 
    
    L_sens = dissipator(gamma_sens, s1, H) +...
         dissipator(gamma_sens, s2, H);

    L = LH + L_atom + L_pump + L_mode + L_sens; 
    
    if j==1
            [A,D] = eigs(L,1,0);    
            %eigs(A,k,sigma) returns k eigenvalues with value=sigma
     else
            [A,D] = eigs(L,1,0,opts);
            %eigs(A,K,sigma,opts) specifies an options structure. 
     end
    
    rho_ss_tmp_func(:,j) = A;
    opts.v0 = A;

end 
   
dY = rho_ss_tmp_func; return

end

dY = L*Y;