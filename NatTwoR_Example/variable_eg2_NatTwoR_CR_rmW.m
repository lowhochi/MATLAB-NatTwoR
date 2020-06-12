% variable_eg2_NatTwoR_CR_rmW.m
muVec = [mu1; mu2; mu3];
conjMuVec = [conj(mu1); conj(mu2); conj(mu3)];
vnormvVec = [v1normv; v2normv; v3normv];
syms x y z real

Mx = [1/(x+y*z), 0, 0;
   0, 1, 0;
   0, 0, 1];

temp0 = (1+u^2)*z + 2*i*u*y;
f = (1-u^2)/(2*(x+y*z))*temp0;
aVs = (1-conj(u)^2)/(2*(x+y*z)*(1+u*conj(u))^2)*temp0;
thetas = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daVs_x = diff(aVs, x);
daVs_y = diff(aVs, y);
daVs_z = diff(aVs, z);
daVs_u = complexdiff3(aVs, u, 0);
daVs_conju = complexdiff3(aVs, u, 1);
dthetas_x = diff(thetas, x);
dthetas_y = diff(thetas, y);
dthetas_z = diff(thetas, z);
daVsVec= [daVs_x, daVs_y, daVs_z]*Mx;
dthetasVec = [dthetas_x, dthetas_y, dthetas_z]*Mx;
% % % % %
daVs_Mu = daVsVec*muVec;
daVs_conjMu = daVsVec*conjMuVec;
daVs_vnormv = daVsVec*vnormvVec;
dthetas_Mu = dthetasVec*muVec;
dthetas_vnormv = dthetasVec*vnormvVec;

daVs_u = complex_simple3(daVs_u, u);
daVs_conju = complex_simple3(daVs_conju, u);
daVs_Mu = complex_simple3(daVs_Mu, u);
daVs_conjMu = complex_simple3(daVs_conjMu, u);
daVs_vnormv = complex_simple3(daVs_vnormv, u);
dthetas_Mu = complex_simple3(dthetas_Mu, u);
dthetas_vnormv = complex_simple3(dthetas_vnormv, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2aVs_uu = complexdiff3(daVs_u, u, 0);
d2aVs_ux = diff(daVs_u, x);
d2aVs_uy = diff(daVs_u, y);
d2aVs_uz = diff(daVs_u, z);
d2aVs_conjux = diff(daVs_conju, x);
d2aVs_conjuy = diff(daVs_conju, y);
d2aVs_conjuz = diff(daVs_conju, z);
d2aVs_uVec = [d2aVs_ux, d2aVs_uy, d2aVs_uz]*Mx;
d2aVs_conjuVec = [d2aVs_conjux, d2aVs_conjuy, d2aVs_conjuz]*Mx;
% % % % %
d2aVs_uMu = d2aVs_uVec*muVec;
d2aVs_uconjMu = d2aVs_uVec*conjMuVec;
d2aVs_uvnormv = d2aVs_uVec*vnormvVec;
d2aVs_conjuMu = d2aVs_conjuVec*muVec;
d2aVs_conjuconjMu = d2aVs_conjuVec*conjMuVec;
d2aVs_conjuvnormv = d2aVs_conjuVec*vnormvVec;

d2aVs_uu = complex_simple3(d2aVs_uu,u);
d2aVs_uMu = complex_simple3(d2aVs_uMu,u);
d2aVs_uconjMu = complex_simple3(d2aVs_uconjMu,u);
d2aVs_uvnormv = complex_simple3(d2aVs_uvnormv,u);
d2aVs_conjuMu = complex_simple3(d2aVs_conjuMu,u);
d2aVs_conjuconjMu = complex_simple3(d2aVs_conjuconjMu,u);
d2aVs_conjuvnormv = complex_simple3(d2aVs_conjuvnormv,u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2aVs_Mux = diff(daVs_Mu, x);
d2aVs_Muy = diff(daVs_Mu, y);
d2aVs_Muz = diff(daVs_Mu, z);
d2aVs_conjMux = diff(daVs_conjMu, x);
d2aVs_conjMuy = diff(daVs_conjMu, y);
d2aVs_conjMuz = diff(daVs_conjMu, z);
d2aVs_vnormvx = diff(daVs_vnormv, x);
d2aVs_vnormvy = diff(daVs_vnormv, y);
d2aVs_vnormvz = diff(daVs_vnormv, z);
d3aVs_uux = diff(d2aVs_uu,x);
d3aVs_uuy = diff(d2aVs_uu,y);
d3aVs_uuz = diff(d2aVs_uu,z);
d2aVs_MuVec = [d2aVs_Mux, d2aVs_Muy, d2aVs_Muz]*Mx;
d2aVs_conjMuVec = [d2aVs_conjMux, d2aVs_conjMuy, d2aVs_conjMuz]*Mx;
d2aVs_vnormvVec = [d2aVs_vnormvx, d2aVs_vnormvy, d2aVs_vnormvz]*Mx;
d3aVs_uuVec = [d3aVs_uux, d3aVs_uuy, d3aVs_uuz]*Mx;
% % % % %
d2aVs_MuMu = d2aVs_MuVec*muVec;
d2aVs_MuconjMu = d2aVs_MuVec*conjMuVec;
d2aVs_Muvnormv = d2aVs_MuVec*vnormvVec;
d2aVs_conjMuconjMu = d2aVs_conjMuVec*conjMuVec;
d2aVs_conjMuvnormv = d2aVs_conjMuVec*vnormvVec;
d2aVs_vnormvvnormv = d2aVs_vnormvVec*vnormvVec;
d3aVs_uuMu = d3aVs_uuVec*muVec;
d3aVs_uuconjMu = d3aVs_uuVec*conjMuVec;
d3aVs_uuvnormv = d3aVs_uuVec*vnormvVec;

d2aVs_MuMu = complex_simple3(d2aVs_MuMu, u);
d2aVs_MuconjMu = complex_simple3(d2aVs_MuconjMu, u);
d2aVs_Muvnormv = complex_simple3(d2aVs_Muvnormv, u);
d2aVs_conjMuconjMu = complex_simple3(d2aVs_conjMuconjMu, u);
d2aVs_conjMuvnormv = complex_simple3(d2aVs_conjMuvnormv, u);
d3aVs_uuMu = complex_simple3(d3aVs_uuMu, u);
d3aVs_uuconjMu = complex_simple3(d3aVs_uuconjMu, u);
d3aVs_uuvnormv = complex_simple3(d3aVs_uuvnormv, u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2thetas_Mux = diff(dthetas_Mu, x);
d2thetas_Muy = diff(dthetas_Mu, y);
d2thetas_Muz = diff(dthetas_Mu, z);
d2thetas_vnormvx = diff(dthetas_vnormv, x);
d2thetas_vnormvy = diff(dthetas_vnormv, y);
d2thetas_vnormvz = diff(dthetas_vnormv, z);
d3aVs_uMux = diff(d2aVs_uMu, x);
d3aVs_uMuy = diff(d2aVs_uMu, y);
d3aVs_uMuz = diff(d2aVs_uMu, z);
d3aVs_uconjMux = diff(d2aVs_uconjMu, x);
d3aVs_uconjMuy = diff(d2aVs_uconjMu, y);
d3aVs_uconjMuz = diff(d2aVs_uconjMu, z);
d3aVs_conjuMux = diff(d2aVs_conjuMu,x);
d3aVs_conjuMuy = diff(d2aVs_conjuMu,y);
d3aVs_conjuMuz = diff(d2aVs_conjuMu,z);
d3aVs_conjuconjMux = diff(d2aVs_conjuconjMu,x);
d3aVs_conjuconjMuy = diff(d2aVs_conjuconjMu,y);
d3aVs_conjuconjMuz = diff(d2aVs_conjuconjMu,z);

d2thetas_MuVec = [d2thetas_Mux, d2thetas_Muy, d2thetas_Muz]*Mx;
d2thetas_vnormvVec = [d2thetas_vnormvx, d2thetas_vnormvy, d2thetas_vnormvz]*Mx;
d3aVs_uMuVec = [d3aVs_uMux, d3aVs_uMuy, d3aVs_uMuz]*Mx;
d3aVs_uconjMuVec = [d3aVs_uconjMux, d3aVs_uconjMuy, d3aVs_uconjMuz]*Mx;
d3aVs_conjuMuVec = [d3aVs_conjuMux, d3aVs_conjuMuy, d3aVs_conjuMuz]*Mx;
d3aVs_conjuconjMuVec = [d3aVs_conjuconjMux, d3aVs_conjuconjMuy,...
    d3aVs_conjuconjMuz]*Mx;
% % % % %
d2thetas_MuMu = d2thetas_MuVec*muVec;
d2thetas_MuconjMu = d2thetas_MuVec*conjMuVec;
d2thetas_Muvnormv = d2thetas_MuVec*vnormvVec;
d2thetas_vnormvvnormv = d2thetas_vnormvVec*vnormvVec;
d3aVs_uMuMu = d3aVs_uMuVec*muVec;
d3aVs_uMuconjMu = d3aVs_uMuVec*conjMuVec;
d3aVs_uMuvnormv = d3aVs_uMuVec*vnormvVec;
d3aVs_uconjMuvnormv = d3aVs_uconjMuVec*vnormvVec;
d3aVs_conjuMuMu = d3aVs_conjuMuVec*muVec;
d3aVs_conjuMuconjMu = d3aVs_conjuMuVec*conjMuVec;
d3aVs_conjuMuvnormv = d3aVs_conjuMuVec*vnormvVec;
d3aVs_conjuconjMuconjMu = d3aVs_conjuconjMuVec*conjMuVec;
d3aVs_conjuconjMuvnormv = d3aVs_conjuconjMuVec*vnormvVec;

d2thetas_MuMu = complex_simple3(d2thetas_MuMu, u);
d2thetas_MuconjMu = complex_simple3(d2thetas_MuconjMu, u);
d2thetas_Muvnormv = complex_simple3(d2thetas_Muvnormv, u);
d2thetas_vnormvvnormv = complex_simple3(d2thetas_vnormvvnormv, u);
d3aVs_uMuMu = complex_simple3(d3aVs_uMuMu, u);
d3aVs_uMuconjMu = complex_simple3(d3aVs_uMuconjMu, u);
d3aVs_uMuvnormv = complex_simple3(d3aVs_uMuvnormv, u);
d3aVs_uconjMuvnormv = complex_simple3(d3aVs_uconjMuvnormv, u);
d3aVs_conjuMuMu = complex_simple3(d3aVs_conjuMuMu, u);
d3aVs_conjuMuconjMu = complex_simple3(d3aVs_conjuMuconjMu, u);
d3aVs_conjuMuvnormv = complex_simple3(d3aVs_conjuMuvnormv, u);
d3aVs_conjuconjMuconjMu = complex_simple3(d3aVs_conjuconjMuconjMu, u);
d3aVs_conjuconjMuvnormv = complex_simple3(d3aVs_conjuconjMuvnormv, u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df_x = diff(f, x);
df_y = diff(f, y);
df_z = diff(f, z);
df_u = complexdiff3(f, u, 0);
dfVec = [df_x, df_y, df_z]*Mx;
% % % % %
df_Mu = dfVec*muVec;
df_conjMu = dfVec*conjMuVec;
df_vnormv = dfVec*vnormvVec;

df_u = complex_simple3(df_u, u);
df_Mu = complex_simple3(df_Mu, u);
df_conjMu = complex_simple3(df_conjMu, u);
df_vnormv = complex_simple3(df_vnormv, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2f_ux = diff(df_u, x);
d2f_uy = diff(df_u, y);
d2f_uz = diff(df_u, z);
d2f_uu = complexdiff3(df_u, u, 0);
d2f_uVec = [d2f_ux, d2f_uy, d2f_uz]*Mx;
% % % % %
d2f_uMu = d2f_uVec*muVec;
d2f_uconjMu = d2f_uVec*conjMuVec;
d2f_uvnormv = d2f_uVec*vnormvVec;

d2f_uu = complex_simple3(d2f_uu, u);
d2f_uMu = complex_simple3(d2f_uMu, u);
d2f_uconjMu = complex_simple3(d2f_uconjMu, u);
d2f_uvnormv = complex_simple3(d2f_uvnormv, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2f_Mux = diff(df_Mu, x);
d2f_Muy = diff(df_Mu, y);
d2f_Muz = diff(df_Mu, z);
d2f_conjMux = diff(df_conjMu, x);
d2f_conjMuy = diff(df_conjMu, y);
d2f_conjMuz = diff(df_conjMu, z);
d2f_vnormvx = diff(df_vnormv, x);
d2f_vnormvy = diff(df_vnormv, y);
d2f_vnormvz = diff(df_vnormv, z);
d3f_uux = diff(d2f_uu, x);
d3f_uuy = diff(d2f_uu, y);
d3f_uuz = diff(d2f_uu, z);
d3f_uuu = complexdiff3(d2f_uu, u, 0);

d2f_MuVec = [d2f_Mux, d2f_Muy, d2f_Muz]*Mx;
d2f_conjMuVec = [d2f_conjMux, d2f_conjMuy, d2f_conjMuz]*Mx;
d2f_vnormvVec = [d2f_vnormvx, d2f_vnormvy, d2f_vnormvz]*Mx;
d3f_uuVec = [d3f_uux, d3f_uuy, d3f_uuz]*Mx;
% % % % %
d2f_MuMu = d2f_MuVec*muVec;
d2f_MuconjMu = d2f_MuVec*conjMuVec;
d2f_Muvnormv = d2f_MuVec*vnormvVec;
d2f_conjMuconjMu = d2f_conjMuVec*conjMuVec;
d2f_conjMuvnormv = d2f_conjMuVec*vnormvVec;
d2f_vnormvvnormv = d2f_vnormvVec*vnormvVec;
d3f_uuMu = d3f_uuVec*muVec;
d3f_uuconjMu = d3f_uuVec*conjMuVec;
d3f_uuvnormv = d3f_uuVec*vnormvVec;

d3f_uuu = complex_simple3(d3f_uuu, u);
d2f_MuMu = complex_simple3(d2f_MuMu, u);
d2f_MuconjMu = complex_simple3(d2f_MuconjMu, u);
d2f_Muvnormv = complex_simple3(d2f_Muvnormv, u);
d2f_conjMuconjMu = complex_simple3(d2f_conjMuconjMu, u);
d2f_conjMuvnormv = complex_simple3(d2f_conjMuvnormv, u);
d2f_vnormvvnormv = complex_simple3(d2f_vnormvvnormv, u);
d3f_uuMu = complex_simple3(d3f_uuMu, u);
d3f_uuconjMu = complex_simple3(d3f_uuconjMu, u);
d3f_uuvnormv = complex_simple3(d3f_uuvnormv, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d3f_uMux = diff(d2f_uMu, x);
d3f_uMuy = diff(d2f_uMu, y);
d3f_uMuz = diff(d2f_uMu, z);
d3f_uconjMux = diff(d2f_uconjMu, x);
d3f_uconjMuy = diff(d2f_uconjMu, y);
d3f_uconjMuz = diff(d2f_uconjMu, z);
d4f_uuMux = diff(d3f_uuMu, x);
d4f_uuMuy = diff(d3f_uuMu, y);
d4f_uuMuz = diff(d3f_uuMu, z);
d4f_uuconjMux = diff(d3f_uuconjMu, x);
d4f_uuconjMuy = diff(d3f_uuconjMu, y);
d4f_uuconjMuz = diff(d3f_uuconjMu, z);
d4f_uuux = diff(d3f_uuu, x);
d4f_uuuy = diff(d3f_uuu, y);
d4f_uuuz = diff(d3f_uuu, z);
d4f_uuuu = complexdiff3(d3f_uuu, u, 0);

d3f_uMuVec = [d3f_uMux, d3f_uMuy, d3f_uMuz]*Mx;
d3f_uconjMuVec = [d3f_uconjMux, d3f_uconjMuy, d3f_uconjMuz]*Mx;
d4f_uuMuVec = [d4f_uuMux, d4f_uuMuy, d4f_uuMuz]*Mx;
d4f_uuconjMuVec = [d4f_uuconjMux, d4f_uuconjMuy, d4f_uuconjMuz]*Mx;
d4f_uuuVec = [d4f_uuux, d4f_uuuy, d4f_uuuz]*Mx;
% % % % %
d3f_uMuMu = d3f_uMuVec*muVec;
d3f_uMuconjMu = d3f_uMuVec*conjMuVec;
d3f_uMuvnormv = d3f_uMuVec*vnormvVec;
d3f_uconjMuconjMu = d3f_uconjMuVec*conjMuVec;
d3f_uconjMuvnormv = d3f_uconjMuVec*vnormvVec;
d4f_uuMuMu = d4f_uuMuVec*muVec;
d4f_uuMuconjMu = d4f_uuMuVec*conjMuVec;
d4f_uuconjMuconjMu = d4f_uuconjMuVec*conjMuVec;
d4f_uuuMu = d4f_uuuVec*muVec;
d4f_uuuconjMu = d4f_uuuVec*conjMuVec;

d4f_uuuu = complex_simple3(d4f_uuuu, u);
d3f_uMuMu = complex_simple3(d3f_uMuMu, u);
d3f_uMuconjMu = complex_simple3(d3f_uMuconjMu, u);
d3f_uMuvnormv = complex_simple3(d3f_uMuvnormv, u);
d3f_uconjMuconjMu = complex_simple3(d3f_uconjMuconjMu, u);
d3f_uconjMuvnormv = complex_simple3(d3f_uconjMuvnormv, u);
d4f_uuMuMu = complex_simple3(d4f_uuMuMu, u);
d4f_uuMuconjMu = complex_simple3(d4f_uuMuconjMu, u);
d4f_uuconjMuconjMu = complex_simple3(d4f_uuconjMuconjMu, u);
d4f_uuuMu = complex_simple3(d4f_uuuMu, u);
d4f_uuuconjMu = complex_simple3(d4f_uuuconjMu, u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%