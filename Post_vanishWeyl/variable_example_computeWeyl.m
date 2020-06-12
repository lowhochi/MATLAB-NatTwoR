% variable_example_computeWeyl.m
daVs_x = diff(aVs, x);
daVs_y = diff(aVs, y);
daVs_z = diff(aVs, z);
daVs_u = complexdiff3(aVs, u, 0);
daVs_conju = complexdiff3(aVs, u, 1);
dthetas_x = diff(thetas, x);
dthetas_y = diff(thetas, y);
dthetas_z = diff(thetas, z);
daVsVec= [daVs_x, daVs_y, daVs_z];
dthetasVec = [dthetas_x, dthetas_y, dthetas_z];
% % % % %
daVs_Mu = daVsVec*muVec;
daVs_conjMu = daVsVec*conjMuVec;
daVs_vnormv = daVsVec*vnormvVec;
dthetas_Mu = dthetasVec*muVec;
dthetas_vnormv = dthetasVec*vnormvVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2aVs_uu = complexdiff3(daVs_u, u, 0);
d2aVs_ux = diff(daVs_u, x);
d2aVs_uy = diff(daVs_u, y);
d2aVs_uz = diff(daVs_u, z);
d2aVs_conjux = diff(daVs_conju, x);
d2aVs_conjuy = diff(daVs_conju, y);
d2aVs_conjuz = diff(daVs_conju, z);
d2aVs_uVec = [d2aVs_ux, d2aVs_uy, d2aVs_uz];
d2aVs_conjuVec = [d2aVs_conjux, d2aVs_conjuy, d2aVs_conjuz];
% % % % %
d2aVs_uMu = d2aVs_uVec*muVec;
d2aVs_uconjMu = d2aVs_uVec*conjMuVec;
d2aVs_uvnormv = d2aVs_uVec*vnormvVec;
d2aVs_conjuMu = d2aVs_conjuVec*muVec;
d2aVs_conjuconjMu = d2aVs_conjuVec*conjMuVec;
d2aVs_conjuvnormv = d2aVs_conjuVec*vnormvVec;
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
d2aVs_MuVec = [d2aVs_Mux, d2aVs_Muy, d2aVs_Muz];
d2aVs_conjMuVec = [d2aVs_conjMux, d2aVs_conjMuy, d2aVs_conjMuz];
d2aVs_vnormvVec = [d2aVs_vnormvx, d2aVs_vnormvy, d2aVs_vnormvz];
d3aVs_uuVec = [d3aVs_uux, d3aVs_uuy, d3aVs_uuz];
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
d2thetas_MuVec = [d2thetas_Mux, d2thetas_Muy, d2thetas_Muz];
d2thetas_vnormvVec = [d2thetas_vnormvx, d2thetas_vnormvy, d2thetas_vnormvz];
d3aVs_uMuVec = [d3aVs_uMux, d3aVs_uMuy, d3aVs_uMuz];
d3aVs_uconjMuVec = [d3aVs_uconjMux, d3aVs_uconjMuy, d3aVs_uconjMuz];
d3aVs_conjuMuVec = [d3aVs_conjuMux, d3aVs_conjuMuy, d3aVs_conjuMuz];
d3aVs_conjuconjMuVec = [d3aVs_conjuconjMux, d3aVs_conjuconjMuy,...
    d3aVs_conjuconjMuz];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df_x = diff(f, x);
df_y = diff(f, y);
df_z = diff(f, z);
df_u = complexdiff3(f, u, 0);
dfVec = [df_x, df_y, df_z];
% % % % %
df_Mu = dfVec*muVec;
df_conjMu = dfVec*conjMuVec;
df_vnormv = dfVec*vnormvVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2f_ux = diff(df_u, x);
d2f_uy = diff(df_u, y);
d2f_uz = diff(df_u, z);
d2f_uu = complexdiff3(df_u, u, 0);
d2f_uVec = [d2f_ux, d2f_uy, d2f_uz];
% % % % %
d2f_uMu = d2f_uVec*muVec;
d2f_uconjMu = d2f_uVec*conjMuVec;
d2f_uvnormv = d2f_uVec*vnormvVec;
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

d2f_MuVec = [d2f_Mux, d2f_Muy, d2f_Muz];
d2f_conjMuVec = [d2f_conjMux, d2f_conjMuy, d2f_conjMuz];
d2f_vnormvVec = [d2f_vnormvx, d2f_vnormvy, d2f_vnormvz];
d3f_uuVec = [d3f_uux, d3f_uuy, d3f_uuz];
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

d3f_uMuVec = [d3f_uMux, d3f_uMuy, d3f_uMuz];
d3f_uconjMuVec = [d3f_uconjMux, d3f_uconjMuy, d3f_uconjMuz];
d4f_uuMuVec = [d4f_uuMux, d4f_uuMuy, d4f_uuMuz];
d4f_uuconjMuVec = [d4f_uuconjMux, d4f_uuconjMuy, d4f_uuconjMuz];
d4f_uuuVec = [d4f_uuux, d4f_uuuy, d4f_uuuz];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSet1 = [aV, daV_u, daV_conju, daV_Mu, daV_conjMu, daV_vnormv,...
    dtheta_Mu, dtheta_vnormv];
variableSet2 = [d2aV_uu, d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv,...
    d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv];
variableSet3 = [d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv, d2aV_conjMuconjMu,...
    d2aV_conjMuvnormv, d2aV_vnormvvnormv, d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv];
variableSet4 = [d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv,...
    d3aV_uconjMuvnormv, d3aV_conjuMuMu, d3aV_conjuMuconjMu,...
    d3aV_conjuMuvnormv, d3aV_conjuconjMuconjMu, d3aV_conjuconjMuvnormv];

variableSetT = [theta, dtheta_Mu, dtheta_vnormv,...
    d2theta_MuMu, d2theta_MuconjMu, d2theta_Muvnormv,...
    d2theta_vnormvvnormv];

variableSetf1 = [w, dw_Mu, dw_conjMu, dw_vnormv, dw_u];
variableSetf2 = [d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu];
variableSetf3 = [d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv,...
    d2w_conjMuconjMu, d2w_conjMuvnormv, d2w_vnormvvnormv,...
    d3w_uuMu, d3w_uuconjMu, d3w_uuu, d3w_uuvnormv];
variableSetf4 = [d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv,...
    d3w_uconjMuconjMu, d3w_uconjMuvnormv,...
    d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuconjMuconjMu,...
    d4w_uuuMu, d4w_uuuconjMu, d4w_uuuu];