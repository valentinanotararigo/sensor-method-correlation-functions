function Ldag = Ldag_gen(w1,w2)

load my_operators.mat;
load parameters.mat;

H_JC = ht * (sigma_p*a + sigma_m*a');
H_int = ht * (epsilon1 * (a*s1'+a'*s1) + epsilon2 * (a*s2'+a'*s2));
H_sens = ht * (w1 * s1' * s1 + w2 * s2' * s2);

H = H_JC + H_int + H_sens;

LH_dag = -commutator(H);

L_atom_dag = dissip_dag(gamma_atom, sigma_m, H);    

L_pump_dag = dissip_dag(gamma_pump, sigma_p, H); 

L_mode_dag = dissip_dag(gamma_mode, a, H); 
    
L_sens_dag = dissip_dag(gamma_sens, s1, H) + dissip_dag(gamma_sens, s2, H);

Ldag = LH_dag + L_atom_dag + L_pump_dag + L_mode_dag + L_sens_dag;    

