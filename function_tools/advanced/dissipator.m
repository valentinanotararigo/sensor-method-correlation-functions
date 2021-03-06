
% Function for dissipator terms to solve ME

function D = dissipator(rate, c, Hamiltonian)

%D = rate * (kron((c').',c) ...
%        - 0.5 * kron((c'*c).',eye(length(Hamiltonian))) ...
%        - 0.5 * kron(eye(length(Hamiltonian)),c'*c));
D = rate * (kron(sparse(c').',sparse(c)) ...
        - 0.5 * kron(sparse(c'*c).',speye(length(Hamiltonian))) ...
        - 0.5 * kron(speye(length(Hamiltonian)),sparse(c'*c)));
end
