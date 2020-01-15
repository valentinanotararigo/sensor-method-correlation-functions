
% function to trace over the total rho and get rho_sensor

function rf = rho_atom_func(r)

evol_bd_op;

r_bis = reshape(r,[size(r,1),sqrt(size(r,2)),sqrt(size(r,2))]);
rf = zeros(size(r,1),length(rho_at_in),length(rho_at_in));

for j = 1:size(r_bis,1)
   
    temporary = squeeze(r_bis(j,:,:));     
    rf(j,:,:) = TrX(temporary,2,[length(rho_at_in),length(rho_mode)]); 
end