
% operators

n=5;   
ntrunc = [1,n-1,1,1];  %the two level atom, the mode, 2 sensors

sigma_p = func_operators(ntrunc,1);
a = func_operators(ntrunc,2);
s1 = func_operators(ntrunc,3);
s2 = func_operators(ntrunc,4);

sigma_m = sigma_p';

save my_operators.mat