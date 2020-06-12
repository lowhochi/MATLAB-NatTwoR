% Compute for the Chrisoffel symbols of the Fefferman metric F0.
function Riem_Gamma=riem_christoffel_Feff_twistor_CR_funW(F0,...
    RVar3, CVar3, dRVar3, dCVar3, MVar3)
% df = df_Feff_twistor_CR_funW(f,RVar3,CVar3,dRVar3,dCVar3,MVar3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differentiate F0(m,n) by (x,y,z,u1,u2,gamma)
% using df_Feff_twistor_CR_funW.m
variable = {'v11','v12','v13','v14','v15','v16';
    'v21','v22','v23','v24','v25','v26';
    'v31','v32','v33','v34','v35','v36';
    'v41','v42','v43','v44','v45','v46';
    'v51','v52','v53','v54','v55','v56';
    'v61','v62','v63','v64','v65','v66'};
for m=1:6
    for n=1:6
        dF0.(variable{m,n})= df_Feff_twistor_CR_funW(F0(m,n),...
            RVar3,CVar3,dRVar3,dCVar3,MVar3);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riem_Gamma(m,n,k) = \Gamma_{mn}^k 
F0_inv = inv(F0);
Riem_Gamma = sym('Riem_Gamma', [6,6,6]);
for m=1:6
    for n=1:6
        for k=1:6
            term = 0;
            for ll=1:6
                term_m = dF0.(variable{ll,n})(m);
                term_n = dF0.(variable{m,ll})(n);
                term_ll = dF0.(variable{m,n})(ll);
                s = term_m + term_n - term_ll;
                term = term + 0.5*F0_inv(k,ll)*s;
                clearvars term_m term_n term_ll s
            end
            term = complex_simple3(term, MVar3);
            Riem_Gamma(m,n,k) = term;
            clear term
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        