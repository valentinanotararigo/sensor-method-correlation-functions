
ht=1;
w_mode = 0;
g=1;

gamma_atom = 0.01*g;
gamma_pump = gamma_atom;
gamma_mode = 0.1*g;
gamma_sens = 2*gamma_mode + gamma_atom;

epsilon1 = 10^-4;
epsilon2 = 10^-4;

R = sqrt(g^2-((gamma_mode-gamma_atom)/4)^2);
R2_minus = sqrt(2*g^2-((gamma_mode-gamma_atom)/4)^2)...
           -sqrt(g^2-((gamma_mode-gamma_atom)/4)^2);
R2_plus = sqrt(2*g^2-((gamma_mode-gamma_atom)/4)^2)...
           + sqrt(g^2-((gamma_mode-gamma_atom)/4)^2);

save parameters.mat