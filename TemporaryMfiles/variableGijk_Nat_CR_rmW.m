syms dG12_3_by1 dG12_3_by2 dG12_3_by3 real
syms dG23_1_by1 dG23_1_by2 dG23_1_by3 real
syms dG31_2_by1 dG31_2_by2 dG31_2_by3 real
syms dG11_2_by1 dG11_2_by2 dG11_2_by3 real
syms dG11_3_by1 dG11_3_by2 dG11_3_by3 real
syms dG22_1_by1 dG22_1_by2 dG22_1_by3 real
syms dG22_3_by1 dG22_3_by2 dG22_3_by3 real
syms dG33_1_by1 dG33_1_by2 dG33_1_by3 real
syms dG33_2_by1 dG33_2_by2 dG33_2_by3 real
v1normv = complex_simple3(v1/norm_of_v,[u]);
v2normv = complex_simple3(v2/norm_of_v,[u]);
v3normv = complex_simple3(v3/norm_of_v,[u]);

dG12_3_MuN = mu1*dG12_3_by1 + mu2*dG12_3_by2 + mu3*dG12_3_by3;
dG12_3_conjMuN=conj(mu1)*dG12_3_by1+conj(mu2)*dG12_3_by2+conj(mu3)*dG12_3_by3;
dG12_3_vnormvN = v1normv*dG12_3_by1+v2normv*dG12_3_by2+v3normv*dG12_3_by3;
dG23_1_MuN = mu1*dG23_1_by1 + mu2*dG23_1_by2 + mu3*dG23_1_by3;
dG23_1_conjMuN=conj(mu1)*dG23_1_by1+conj(mu2)*dG23_1_by2+conj(mu3)*dG23_1_by3;
dG23_1_vnormvN = v1normv*dG23_1_by1+v2normv*dG23_1_by2+v3normv*dG23_1_by3;
dG31_2_MuN = mu1*dG31_2_by1 + mu2*dG31_2_by2 + mu3*dG31_2_by3;
dG31_2_conjMuN=conj(mu1)*dG31_2_by1+conj(mu2)*dG31_2_by2+conj(mu3)*dG31_2_by3;
dG31_2_vnormvN = v1normv*dG31_2_by1+v2normv*dG31_2_by2+v3normv*dG31_2_by3;
dG11_2_MuN = mu1*dG11_2_by1 + mu2*dG11_2_by2 + mu3*dG11_2_by3;
dG11_2_conjMuN=conj(mu1)*dG11_2_by1+conj(mu2)*dG11_2_by2+conj(mu3)*dG11_2_by3;
dG11_2_vnormvN = v1normv*dG11_2_by1+v2normv*dG11_2_by2+v3normv*dG11_2_by3;
dG11_3_MuN = mu1*dG11_3_by1 + mu2*dG11_3_by2 + mu3*dG11_3_by3;
dG11_3_conjMuN=conj(mu1)*dG11_3_by1+conj(mu2)*dG11_3_by2+conj(mu3)*dG11_3_by3;
dG11_3_vnormvN = v1normv*dG11_3_by1+v2normv*dG11_3_by2+v3normv*dG11_3_by3;
dG22_1_MuN = mu1*dG22_1_by1 + mu2*dG22_1_by2 + mu3*dG22_1_by3;
dG22_1_conjMuN=conj(mu1)*dG22_1_by1+conj(mu2)*dG22_1_by2+conj(mu3)*dG22_1_by3;
dG22_1_vnormvN = v1normv*dG22_1_by1+v2normv*dG22_1_by2+v3normv*dG22_1_by3;
dG22_3_MuN = mu1*dG22_3_by1 + mu2*dG22_3_by2 + mu3*dG22_3_by3;
dG22_3_conjMuN=conj(mu1)*dG22_3_by1+conj(mu2)*dG22_3_by2+conj(mu3)*dG22_3_by3;
dG22_3_vnormvN = v1normv*dG22_3_by1+v2normv*dG22_3_by2+v3normv*dG22_3_by3;
dG33_1_MuN = mu1*dG33_1_by1 + mu2*dG33_1_by2 + mu3*dG33_1_by3;
dG33_1_conjMuN=conj(mu1)*dG33_1_by1+conj(mu2)*dG33_1_by2+conj(mu3)*dG33_1_by3;
dG33_1_vnormvN = v1normv*dG33_1_by1+v2normv*dG33_1_by2+v3normv*dG33_1_by3;
dG33_2_MuN = mu1*dG33_2_by1 + mu2*dG33_2_by2 + mu3*dG33_2_by3;
dG33_2_conjMuN=conj(mu1)*dG33_2_by1+conj(mu2)*dG33_2_by2+conj(mu3)*dG33_2_by3;
dG33_2_vnormvN = v1normv*dG33_2_by1+v2normv*dG33_2_by2+v3normv*dG33_2_by3;

% dG12_3Row = [dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv]
dG12_3N = [dG12_3_MuN, dG12_3_conjMuN, dG12_3_vnormvN];
dG23_1N = [dG23_1_MuN, dG23_1_conjMuN, dG23_1_vnormvN];
dG31_2N = [dG31_2_MuN, dG31_2_conjMuN, dG31_2_vnormvN];
dG11_2N = [dG11_2_MuN, dG11_2_conjMuN, dG11_2_vnormvN];
dG11_3N = [dG11_3_MuN, dG11_3_conjMuN, dG11_3_vnormvN];
dG22_1N = [dG22_1_MuN, dG22_1_conjMuN, dG22_1_vnormvN];
dG22_3N = [dG22_3_MuN, dG22_3_conjMuN, dG22_3_vnormvN];
dG33_1N = [dG33_1_MuN, dG33_1_conjMuN, dG33_1_vnormvN];
dG33_2N = [dG33_2_MuN, dG33_2_conjMuN, dG33_2_vnormvN];