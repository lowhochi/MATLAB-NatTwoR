% flatM2_rewrite_NatTwoR_CR_rmW.m
load('DataFlatM1_NatTwoR_CR_rmW.mat');

d2w_vnormvMu = subs(d2w_vnormvMu, [aV,bV,T4], zeros(1,3));
d2w_vnormvu = d2w_uvnormv + 1/y^2*conj(GH11_2);
d2rho_vnormvMu = subs(d2rho_vnormvMu, [aV,bV,T4], zeros(1,3));

syms dphiW_Mu dphiW_conjMu dphiW_vnormv dphiW_u
dphiW_conju = -6/y^2*conj(GH11_1)+12*conj(u)/y^3*w;

syms dGH11_2_Mu dGH11_2_conjMu dGH11_2_vnormv dGH11_2_conju
dGH11_2_u = -2*conj(u)/y*GH11_1  -4*u*conj(u)/y^2*conj(w)...
    -2*conj(dw_vnormv);

dGT01_2_Mu = w/y^2*GH11_2 -conj(w)/y^2*dw_Mu;
dGT01_2_conjMu = conj(dGT01_2_Mu);
dGT01_2_vnormv = -1/y^2*(conj(w)*dw_vnormv + w*conj(dw_vnormv));
dGT01_2_u = -1/y^2*conj(w)*conj(GH11_1);
dGT01_2_conju = conj(dGT01_2_u);

d2GH11_1_Muu = 2/y^2*GH11_2 +2*conj(u)/y*dGH11_1_Mu...
    +2*dGH11_1_vnormv;
d2GH11_1_conjMuu = -2/y^2*conj(dw_Mu);
d2GH11_1_vnormvu = -2/y^2*conj(dw_vnormv) -1/y^2*dGH11_1_conjMu;

dGH11_1_conjuZ = conj(phiW) +4*u/y*GH11_1 -2*u^2/y^2*conj(w);
d2GH11_1_conjuMuZ = conj(dphiW_conjMu) +4*u/y*dGH11_1_Mu...
    +2*u^2/y^2*GH11_2;
d2GH11_1_conjuconjMuZ = conj(dphiW_Mu) +4*u/y*dGH11_1_conjMu...
    -2*u^2/y^2*conj(dw_Mu);
d2GH11_1_conjuvnormvZ = conj(dphiW_vnormv) +4*u/y*dGH11_1_vnormv...
    -2*u^2/y^2*conj(dw_vnormv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarFlatM2 = [GH11_1, GH11_2, GT01_2, dGH11_1_Mu, dGH11_1_conjMu,...
    dGH11_1_vnormv, drho_Mu, drho_conjMu, drho_conju,...
    drho_u, drho_vnormv, dw_vnormv, phiW, rho, u, w];

% derivativeDict in [D_Mu; D_conjMu; D_vnormv; d/du; d/dconj(u)];
% Redefine 2nd derivatives for the flat case.
derivativeDict.u =[0; 0; 0; 1; 0];
derivativeDict.w =[dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0];
derivativeDict.dw_vnormv = [d2w_vnormvMu; d2w_vnormvconjMu;
    d2w_vnormvvnormv; d2w_vnormvu; d2w_vnormvconju];

derivativeDict.rho=[drho_Mu; conj(drho_Mu); drho_vnormv; drho_u; conj(drho_u)];
derivativeDict.drho_Mu =[d2rho_MuMu; d2rho_MuconjMu; d2rho_Muvnormv;
    d2rho_Muu; d2rho_Muconju];
derivativeDict.drho_conjMu =[conj(d2rho_MuconjMu); conj(d2rho_MuMu);
    conj(d2rho_Muvnormv); conj(d2rho_Muconju); conj(d2rho_Muu)];
derivativeDict.drho_vnormv=[d2rho_vnormvMu; conj(d2rho_vnormvMu);
    d2rho_vnormvvnormv; d2rho_vnormvu; conj(d2rho_vnormvu)];
derivativeDict.drho_u =[d2rho_uMu; d2rho_uconjMu; d2rho_uvnormv; 
    d2rho_uu; d2rho_uconju];
derivativeDict.drho_conju =[conj(d2rho_uconjMu); conj(d2rho_uMu); conj(d2rho_uvnormv); 
    d2rho_uconju; conj(d2rho_uu)];

derivativeDict.phiW = [dphiW_Mu; dphiW_conjMu; dphiW_vnormv;
    dphiW_u; dphiW_conju];
derivativeDict.GH11_1 = [dGH11_1_Mu; dGH11_1_conjMu; dGH11_1_vnormv;
    dGH11_1_u; dGH11_1_conju];
derivativeDict.GH11_2 = [dGH11_2_Mu; dGH11_2_conjMu; dGH11_2_vnormv;
    dGH11_2_u; dGH11_2_conju];
derivativeDict.GT01_2 = [dGT01_2_Mu; dGT01_2_conjMu; dGT01_2_vnormv;
    dGT01_2_u; dGT01_2_conju];
derivativeDict.dGH11_1_Mu = [d2GH11_1_MuMu; d2GH11_1_MuconjMu; d2GH11_1_Muvnormv;
    d2GH11_1_Muu; d2GH11_1_Muconju];
derivativeDict.dGH11_1_conjMu = [d2GH11_1_conjMuMu; d2GH11_1_conjMuconjMu; 
    d2GH11_1_conjMuvnormv; d2GH11_1_conjMuu; d2GH11_1_conjMuconju];
derivativeDict.dGH11_1_vnormv = [d2GH11_1_vnormvMu; d2GH11_1_vnormvconjMu; 
    d2GH11_1_vnormvvnormv; d2GH11_1_vnormvu; d2GH11_1_vnormvconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write down a 1-form in the format,
% [\theta^1, \theta^\bar1, \theta^2, \theta^\bar2, \alpha, d\gamma];
OmegaHat = cell(6,6); %OmegaHat{m,n} = \hat\omega_m^n

for m=1:6
    for n=1:6        
        temp1 = RGammaFlat(1,m,n);
        temp2 = RGammaFlat(2,m,n);    
        temp3 = RGammaFlat(3,m,n);
        temp4 = RGammaFlat(4,m,n); 
        temp5 = RGammaFlat(5,m,n);
        temp6 = RGammaFlat(6,m,n);

        OmegaHat{m,n} = [temp1, temp2, temp3, temp4, temp5, temp6];
    end
end
clearvars m n temp1 temp2 temp3 temp4 temp5 tem6

dtheta1 = 1/2*[0, 0, 0, 2*u/y, 0, 0;
    0, 0, 0, 0, w/y^2, 0;
    0, 0, 0, 0, 1/y^2, 0;
    -2*u/y, 0, 0, 0, 0, 0;
    0, -w/y^2, -1/y^2, 0, 0, 0;
    0, 0, 0, 0, 0, 0];

dtheta1bar = 1/2*[0, 0, 0, 0, conj(w)/y^2, 0;
    0, 0, 2*conj(u)/y, 0, 0, 0;
    0, -2*conj(u)/y, 0, 0, 0, 0;
    0, 0, 0, 0, 1/y^2, 0;
    -conj(w)/y^2, 0, 0, -1/y^2, 0, 0;
    0, 0, 0, 0, 0, 0];

dtheta2 = 1/2*[0, conj(GH11_2), 0, 0, GT01_2, 0;
    -conj(GH11_2), 0, conj(GH11_1), 0, dw_vnormv, 0;
    0, -conj(GH11_1), 0, 0, 0, 0;
    0, 0, 0, 0, -w/y^2, 0;
    -GT01_2, -dw_vnormv, 0, w/y^2, 0, 0;
    0, 0, 0, 0, 0, 0];

dtheta2bar = 1/2*[0, -GH11_2, 0, GH11_1, conj(dw_vnormv), 0;
    GH11_2, 0, 0, 0, GT01_2, 0;
    0, 0, 0, 0, -conj(w)/y^2, 0;
    -GH11_1, 0, 0, 0, 0, 0;
    -conj(dw_vnormv), -GT01_2, conj(w)/y^2, 0 ,0, 0;
    0, 0, 0, 0, 0, 0];

dalpha = [0, 0, 0, 1, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, -1, 0, 0, 0, 0;
    -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0];
dCoframe = {dtheta1, dtheta1bar, dtheta2, dtheta2bar, dalpha, zeros(6,6)};

dOmegaHat = cell(6,6);
for m=1:6
    for n=1:6
        dOmegaHat{m,n} = sym('dOmegaHat',[1,6]);
    end
end

MVarFlat2 = [];
for m=1:6
    for n=1:6
        temp1 = OmegaHat{m,n};
        dTemp = df_flatM3_CR_rmW(temp1, dCoframe, CVarFlatM2, derivativeDict, gamma);
        MVarFlat2 = union(MVarFlat2, symvar(dTemp));
        for j=1:6
            for k=1:6
            temp2 = dTemp(j,k);
            temp2 = subs(temp2, d2GH11_1_conjuMu, d2GH11_1_conjuMuZ);
            temp2 = subs(temp2, d2GH11_1_conjuconjMu, d2GH11_1_conjuconjMuZ);
            temp2 = subs(temp2, d2GH11_1_conjuvnormv, d2GH11_1_conjuvnormvZ);   
            temp2 = subs(temp2, dGH11_1_conju, dGH11_1_conjuZ);
            
            temp2 = subs(temp2, [d3w_uuu, d3w_uuconjMu, d3w_uuMu, d3w_uuvnormv],...
                [d3w_uuuZ, d3w_uuconjMuZ, d3w_uuMuZ, d3w_uuvnormvZ]);
            temp2 = subs(temp2, [d2w_uMu, d2w_uconjMu, d2w_uvnormv],...
                [d2w_uMuZ, d2w_uconjMuZ, d2w_uvnormvZ]);
            temp2 = subs(temp2, [d2w_uu, dw_conjMu], [d2w_uuZ, dw_conjMuZ]);
    
            temp2 = subs(temp2, conj(phiW), phiW+i*rho);
            temp2 = subs(temp2, dw_u, conj(GH11_1)+2*conj(u)/y*w);
            temp2 = subs(temp2, w*conj(w), -y^2*GT01_2);
            dOmegaHat{m,n}(j,k) = temp2;
            end
        end
    end
end

MVarFlat2 = union(MVarFlat2, [Y]);
for m=1:6
    for n=1:6
        for j=1:6
            for k=1:6
            temp = dOmegaHat{m,n}(j,k);
            temp = subs(temp, u*conj(u), Y-1);
            temp = complex_simple3(temp, MVarFlat2);     
            temp = subs(temp, u*conj(u), Y-1);
            dOmegaHat{m,n}(j,k) = complex_simple3(temp, MVarFlat2);
            end
        end
    end
end
clearvars temp temp1 temp2 m n j dTemp

save('DataTemp_FlatM2_NatTwoR_CR_rmW.mat');
