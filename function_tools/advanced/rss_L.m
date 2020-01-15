
% steady state from diagonalization of the Liouvillian

function rss_norm = rss_L(L) 

[A_ss,DD_ss] = eigs(L,1,0);    
opts.v0 = A_ss; %opts.v0=starting vector
[A_ss,DD_ss] = eigs(L,1,0,opts);
%eigs(A,k,sigma) returns k eigenvalues with value=sigma

rss_tmp = A_ss;

%opts.v0 = A_ss;
rss_matr = reshape(rss_tmp,sqrt(length(rss_tmp)),sqrt(length(rss_tmp)));
rss_norm = rss_matr / trace(rss_matr);   % matrix
