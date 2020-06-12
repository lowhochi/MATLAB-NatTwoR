function ricci_curv = ricci_twistor_CR_model(R, h, M)
% ric(i,j) = ric(Xi, conj_Xj)
ric = sym('ric_%d_%d', [2,2]);
rho = 0;
for m=1:2
    for n=1:2
        temp = 0;
        for k=1:2
            temp = temp + R(m,k,k,n);
        end
        ric(m,n) = temp;
        clear temp
    end
end
%
h_inv = inv(h);
for m = 1:2
    for n = 1:2
        rho = rho + h_inv(n,m)*ric(m,n); %Change in rho here.
    end
end
ricci_curv = {ric, rho};

            
           