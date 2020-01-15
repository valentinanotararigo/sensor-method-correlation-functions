
% operators for a three level system, the mode, sensor

n=7;   % to calculate sigma_z I need a_tmp
ntrunc = [1,n-1,1,1]; 
    %the two level atom, the mode, 2 senosors

sigma_m = func_operators(ntrunc,1);
a = func_operators(ntrunc,2);
s1 = func_operators(ntrunc,3);
s2 = func_operators(ntrunc,4);

sigma_p = sigma_m';

%save my_operators.mat