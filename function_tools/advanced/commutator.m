
% Function for the coherent term to solve ME

function L_Hamilt = commutator(Hamilt)

L_Hamilt = 1i * (kron(speye(length(Hamilt)),Hamilt) - ...
                  kron(Hamilt.',speye(length(Hamilt))));

end