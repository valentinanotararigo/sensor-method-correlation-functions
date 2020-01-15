% Function to get the Liouvillian and to solve the Linblad master equation 
% of a two level system coupled to a resonant mode, considering also the 
% interaction with 2 sensors.

function dY = func_2sensors_jk(T,Y)

persistent rho_ss_tmp_func
%persistent L

if nargin==0

parameters;
my_operators;

H_JC = ht * (sigma_p*a + sigma_m*a');
H_int = ht * (epsilon1 * (a*s1'+a'*s1) + epsilon2 * (a*s2'+a'*s2));
H_inv = H_JC + H_int;

L_atom = dissipator(gamma_atom, sigma_m, H_inv);    
L_pump = dissipator(gamma_pump, sigma_p, H_inv);
L_inv = L_atom + L_pump;

for j = 1:length(gamma_mode)
    
    L_mode = dissipator(gamma_mode(j), a, H_inv);

    for k = 1:nloop
    
        H_sens = ht * (w_sens(k) * s1' * s1 + w_2 * s2' * s2);
        H = H_inv + H_sens;
    
        LH = commutator(H);
        L_sens = dissipator(gamma_sens1(j), s1, H_inv) ...
                 + dissipator(gamma_sens2(j), s2, H_inv);
             
        L = LH + L_inv + L_mode + L_sens;
        %L(:,:,j,k) = LH + L_inv + L_mode + L_sens;
        
        if k==1
            [A,D] = eigs(L,1,0);    
            %eigs(A,k,sigma) returns k eigenvalues with value=sigma
        else
            [A,D] = eigs(L,1,0,opts);
            %eigs(A,K,sigma,opts) specifies an options structure. 
        end
    
        rho_ss_tmp_func(:,j,k) = A;
        opts.v0 = A;           %opts.v0=starting vector

    end
end

dY = rho_ss_tmp_func; return
%dY = L; return

end

dY = L*Y;
