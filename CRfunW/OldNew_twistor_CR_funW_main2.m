% New_twistor_CR_funW_main2.m
load('New_twistor_CR_funW_data_Jan17.mat');
assumeAlso([x y z u1 u2 gamma],'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d3w_xxu d3w_xyu d3w_xzu d3w_yyu d3w_yzu d3w_zzu
syms d3w_xuu d3w_yuu d3w_zuu d3w_uuu
syms d4w_xxuu d4w_yyuu d4w_zzuu d4w_xyuu d4w_xzuu d4w_yzuu
syms d4w_xuuu d4w_yuuu d4w_zuuu d4w_uuuu
RVar3 = [x,y,z,u1,u2,gamma];
dRVar3 = eye(6);
CVar3 = [w, dw_u, d2w_uu];
dCVar3 = [dw_x, d2w_xu, d3w_xuu;
    dw_y, d2w_yu, d3w_yuu;
    dw_z, d2w_zu, d3w_zuu;
    dw_u, d2w_uu, d3w_uuu;
    i*dw_u, i*d2w_uu, i*d3w_uuu;
    0, 0, 0];
MVar3 = [u, w, dw_x, dw_y, dw_z, dw_u, d2w_xu, d2w_yu,...
    d2w_zu, d2w_uu, d3w_xuu, d3w_yuu, d3w_zuu, d3w_uuu];
%
Riem_Gamma = riem_christoffel_Feff_twistor_CR_funW(F0,...
    RVar3, CVar3, dRVar3, dCVar3, MVar3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('New_twistor_CR_funW_RiemGamma_data_Jan17.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
