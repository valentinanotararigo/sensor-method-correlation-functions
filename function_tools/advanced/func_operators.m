
function y = func_operators(nt,nn) 

l = 1:nt(nn)+1;

operator = sparse(l,l+1,sqrt(l));
operator(:,length(operator))=[];         %to eliminate the last column 

y = 1;

for npos = 1:length(nt)
    
        if npos == nn;
            v{npos} = operator;
        end
        
        if npos ~= nn;
            v{npos} = eye(nt(npos)+1);
        end
        
        y = kron(y,v{npos});
end

end





