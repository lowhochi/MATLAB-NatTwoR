% variable_Conformal_main1.m
syms u w f %f is real-valued
syms dw_Mu dw_conjMu dw_vnormv dw_u
syms df_Mu df_vnormv df_u
df_conjMu = conj(df_Mu);
df_conju = conj(df_u);
mu1 =u^2-1; 
mu2= 2*u; 
mu3= i*(u^2+1);
v1 = i*(mu2*conj(mu3)-mu3*conj(mu2));
v2 = i*(mu3*conj(mu1)-mu1*conj(mu3));
v3 = i*(mu1*conj(mu2)-mu2*conj(mu1));
v1 = complex_simple3(v1,[u]);
v2 = complex_simple3(v2,[u]);
v3 = complex_simple3(v3,[u]);
norm_of_v = sqrt(v1*v1+v2*v2+v3*v3);
norm_of_v = complex_simple3(norm_of_v,[u]);
v1normv = v1/norm_of_v; %T1
v2normv = v2/norm_of_v; %T2
v3normv = v3/norm_of_v; %T3
v1normv = complex_simple3(v1normv, [u,w]);
v2normv = complex_simple3(v2normv, [u,w]);
v3normv = complex_simple3(v3normv, [u,w]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lieMu = cell(5,5); % [Mu, conjMu, vnormv];
lie_Mu_conjMu = [0; 0; 0];
lie_Mu_vnormv = [0; 0; 0];
lie_Mu_ddu = [-2*conj(u)/(1+u*conj(u)); 0; -2];
lie_conjMu_vnormv = [0; 0; 0];
lie_conjMu_ddconju = [0; -2*u/(1+u*conj(u)); -2];
lie_vnormv_ddu = [0; 1/(1+u*conj(u))^2; 0];
lie_vnormv_ddconju = [1/(1+u*conj(u))^2; 0; 0];
lieMu{1,1} = zeros(3,1);
lieMu{1,2} = lie_Mu_conjMu;
lieMu{1,3} = lie_Mu_vnormv;
lieMu{1,4} = lie_Mu_ddu;
lieMu{1,5} = zeros(3,1);
lieMu{2,1} = -lie_Mu_conjMu;
lieMu{2,2} = zeros(3,1);
lieMu{2,3} = lie_conjMu_vnormv;
lieMu{2,4} = zeros(3,1);
lieMu{2,5} = lie_conjMu_ddconju;
lieMu{3,1} = -lie_Mu_vnormv;
lieMu{3,2} = -lie_conjMu_vnormv;
lieMu{3,3} = zeros(3,1);
lieMu{3,4} = lie_vnormv_ddu;
lieMu{3,5} = lie_vnormv_ddconju;
lieMu{4,1} = -lie_Mu_ddu;
lieMu{4,2} = zeros(3,1);
lieMu{4,3} = -lie_vnormv_ddu;
lieMu{4,4} = zeros(3,1);
lieMu{4,5} = zeros(3,1);
lieMu{5,1} = zeros(3,1);
lieMu{5,2} = -lie_conjMu_ddconju;
lieMu{5,3} = -lie_vnormv_ddconju;
lieMu{5,4} = zeros(3,1);
lieMu{5,5} = zeros(3,1);

syms d2f_MuMu d2f_MuconjMu d2f_Muvnormv d2f_vnormvvnormv 
d2f_conjMuconjMu = conj(d2f_MuMu);
d2f_conjMuvnormv = conj(d2f_Muvnormv);
syms d2f_uu d2f_uconju 
d2f_conjuconju = conj(d2f_uu);
syms d2f_uMu d2f_uconjMu d2f_uvnormv 
d2f_conjuMu = conj(d2f_uconjMu);
d2f_conjuconjMu = conj(d2f_uMu);
d2f_conjuvnormv = conj(d2f_uvnormv);

dfRow = [df_Mu, df_conjMu, df_vnormv];
d2f_conjMuMu = d2f_MuconjMu + dfRow*lieMu{1,2};
d2f_vnormvMu = d2f_Muvnormv + dfRow*lieMu{1,3}; 
d2f_vnormvconjMu = d2f_conjMuvnormv + dfRow*lieMu{2,3};
d2f_Muu = d2f_uMu + dfRow*lieMu{4,1}; 
d2f_conjMuu = d2f_uconjMu + dfRow*lieMu{4,2};
d2f_vnormvu = d2f_uvnormv + dfRow*lieMu{4,3};
d2f_Muconju = d2f_conjuMu + dfRow*lieMu{5,1};
d2f_conjMuconju = d2f_conjuconjMu + dfRow*lieMu{5,2};
d2f_vnormvconju = d2f_conjuvnormv + dfRow*lieMu{5,3};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVar1 = [u, w, f, df_Mu, df_vnormv, df_u];

MVar1 = [u, w, f, dw_Mu, dw_conjMu, dw_vnormv, dw_u,...
    df_Mu, df_vnormv, df_u, d2f_MuMu, d2f_MuconjMu, d2f_Muvnormv,...
    d2f_vnormvvnormv, d2f_uu, d2f_uconju, d2f_uMu, d2f_uconjMu, d2f_uvnormv]; 


derivativeDict.u = [0; 0; 0; 1; 0];
derivativeDict.w = [dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0];
derivativeDict.f = [df_Mu; df_conjMu; df_vnormv; df_u; df_conju];
derivativeDict.df_Mu = [d2f_MuMu; d2f_MuconjMu; d2f_Muvnormv; d2f_Muu; d2f_Muconju];
derivativeDict.df_vnormv = [d2f_vnormvMu; d2f_vnormvconjMu; 
    d2f_vnormvvnormv; d2f_vnormvu; d2f_vnormvconju];
derivativeDict.df_u = [d2f_uMu; d2f_uconjMu; d2f_uvnormv; d2f_uu; d2f_uconju];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
