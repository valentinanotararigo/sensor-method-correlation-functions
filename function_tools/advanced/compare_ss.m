
% comparison between the ss calculated from the Liouvillian and the
% dynamics for long time

function compare_ss(rL,rdyn)

disp('configuration');
rL_vec = reshape(rL,[],1); 
rdyn_vec = reshape(rdyn,[],1);
for k = 1:length(rL_vec)
    diff(k) = rL_vec(k) - rdyn_vec(k);
    if diff(k) > 10^-5
        disp(['diff [',num2str(k),'] = ',num2str(diff(k))]);
    end
end
