function L = L_gen(w1,w2)

load my_operators.mat;
load parameters.mat;

H_JC = ht * (sigma_p*a + sigma_m*a');
H_int = ht * (epsilon1 * (a*s1'+a'*s1) + epsilon2 * (a*s2'+a'*s2));
H_sens = ht * (w1 * s1' * s1 + w2 * s2' * s2);
    
H = H_JC + H_int + H_sens;
    
LH = commutator(H);
L_atom = dissipator(gamma_atom, sigma_m, H);    
L_pump = dissipator(gamma_pump, sigma_p, H); 
L_mode = dissipator(gamma_mode, a, H); 
L_sens = dissipator(gamma_sens, s1, H) + dissipator(gamma_sens, s2, H);

L = LH + L_atom + L_pump + L_mode + L_sens; 