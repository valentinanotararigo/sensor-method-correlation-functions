
% Parameters w_mode = 0;  

w_mode = 0; 
ht=1;
g = 1;

gamma_atom = 0.01*g;
gamma_pump = gamma_atom;

epsilon1 = 0.00001;
epsilon2 = 0.00001;

gamma_mode = [0.01; 0.1; 0.5]*g;

gamma_sens1 = [2*gamma_mode + gamma_atom,...
       (2*gamma_mode + gamma_atom)/2,(2*gamma_mode + gamma_atom)/2,];
gamma_sens2 = gamma_sens1;

nloop = 200;
diff = linspace(-4.5,4.5,nloop);   %diff=(w_sens-w_mode)/g
w_sens = g.*diff + w_mode;

for j = 1:length(gamma_mode)
    
    w_2 = sqrt(g^2-((gamma_mode(j)-gamma_atom)/4)^2);
end

%save parameters.mat
