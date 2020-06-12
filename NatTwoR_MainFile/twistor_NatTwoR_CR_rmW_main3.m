% twistor_NatTwoR_CR_rmW_main3.m
load('DataMain2_NatTwoR_CR_rmW.mat');
assumeAlso(gamma, 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable_NatTwoR_Main3_CR_rmW

dRGamma = cell(6,6,6);
for j=1:6
    for k=1:6
        for ll=1:6
            dRGamma{j,k,ll} = df_NatTwo_MuGamma_CR_rmW(RGamma(j,k,ll),...
                CVarMain3, dCVarMain3, gamma);
            % df = [df/dMu, df/dconjMu, df/dvnormv, df/du, df/dconju, df/dgamma].
        end
    end
end
% RiemCurv0(i,j,k,l) = R_{ijk}^l
% RiemCurv(m,n,k,ll) = R_{mnkll}
RiemCurv0 = sym('RiemCurv0',[6 6 6 6]); 
RiemCurv = sym('RiemCurv',[6 6 6 6]);
symSetRCurv = [];
countIndex1296 = zeros(1296,4);
RowNumber = 1;
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                countIndex1296(RowNumber,:) = [m,n,k,ll];
                RowNumber = RowNumber + 1;
            end
        end
    end
end

for j=1:1296
    m = countIndex1296(j,1);
    n = countIndex1296(j,2);
    k = countIndex1296(j,3);
    ll = countIndex1296(j,4);   
    temp1 = dRGamma{n,k,ll}*Uvector(:,m)-dRGamma{m,k,ll}*Uvector(:,n);
    temp2 = 0;
    for p=1:6
        temp2 = temp2 + RGamma(n,k,p)*RGamma(m,p,ll)...
            -RGamma(m,k,p)*RGamma(n,p,ll)...
            -lieTwo{m,n}(p)*RGamma(p,k,ll);
    end
    RiemCurv0(m,n,k,ll) = temp1 + temp2;  
end
% R_{mnkll} = R_{mnk}^p*F0(p,ll)
for j=1:1296
    m = countIndex1296(j,1);
    n = countIndex1296(j,2);
    k = countIndex1296(j,3);
    ll = countIndex1296(j,4); 
    temp3 = 0;
    for p=1:6
        temp3 = temp3 + F0(p,ll)*RiemCurv0(m,n,k,p);
    end
    RiemCurv(m,n,k,ll) = temp3; 
    tempSet = symvar(temp3);
    symSetRCurv = union(symSetRCurv, tempSet);
end

save('DataTempMain3_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load('DataTempMain3_NatTwoR_CR_rmW.mat');
assumeAlso(gamma, 'real');
RiemCurv_rawVersion = RiemCurv; % unsimplified version of RiemCurv

% Use symmetry to simplify RiemCurv
for j=1:1296
    m = countIndex1296(j,1);
    n = countIndex1296(j,2);
    k = countIndex1296(j,3);
    ll = countIndex1296(j,4); 
    if (m==n)||(k==ll)
        RiemCurv(m,n,k,ll) = 0;
        continue
    end
    
    if (m>n) % Riem(n,m,k,ll) has already been simplified
        RiemCurv(m,n,k,ll) = -RiemCurv(n,m,k,ll);
        continue
    end
    
    if (k>ll) %k>ll with m<n
        RiemCurv(m,n,k,ll) = -RiemCurv(m,n,ll,k);
        continue
    end  
    temp4 = RiemCurv(m,n,k,ll);
    temp4 = complex_simple3(temp4, symSetRCurv);
    RiemCurv(m,n,k,ll) = temp4;
end

clearvars temp1 temp2 temp3 temp4 tempSet
clearvars m n k ll j p RowNumber
clearvars dwRow dh11Row dT4Row daMRow daVRow drhoRow d2wRow d2T4Row
clearvars d2w d2h11 d2T4 d2aM d2aV d2rho d3w_u d3T4_u

list_of_variables_Main3 = who;
save('DataMain3_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%