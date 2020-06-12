% New_twistor_CR_funW_main3.m
load('New_twistor_CR_funW_RiemGamma_data_Jan17.mat');
assumeAlso([x y z u1 u2 gamma],'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RVar4 = [x,y,z,u1,u2,gamma];
dRVar4 = eye(6);
CVar4 = [w, dw_x, dw_y, dw_z, dw_u, d2w_xu, d2w_yu, d2w_zu, d2w_uu,...
    d3w_xuu, d3w_yuu, d3w_zuu, d3w_uuu];
dCVar4 = sym('dCVar4',[6,13]);
dCVar4(:,1)=[dw_x; dw_y; dw_z; dw_u; i*dw_u; 0];
dCVar4(:,2)=[d2w_xx; d2w_xy; d2w_xz; d2w_xu; i*d2w_xu; 0];
dCVar4(:,3)=[d2w_xy; d2w_yy; d2w_yz; d2w_yu; i*d2w_yu; 0];
dCVar4(:,4)=[d2w_xz; d2w_yz; d2w_zz; d2w_zu; i*d2w_zu; 0];
dCVar4(:,5)=[d2w_xu; d2w_yu; d2w_zu; d2w_uu; i*d2w_uu; 0];
dCVar4(:,6)=[d3w_xxu; d3w_xyu; d3w_xzu; d3w_xuu; i*d3w_xuu; 0];
dCVar4(:,7)=[d3w_xyu; d3w_yyu; d3w_yzu; d3w_yuu; i*d3w_yuu; 0];
dCVar4(:,8)=[d3w_xzu; d3w_yzu; d3w_zzu; d3w_zuu; i*d3w_zuu; 0];
dCVar4(:,9)=[d3w_xuu; d3w_yuu; d3w_zuu; d3w_uuu; i*d3w_uuu; 0];
dCVar4(:,10)=[d4w_xxuu; d4w_xyuu; d4w_xzuu; d4w_xuuu; i*d4w_xuuu; 0];
dCVar4(:,11)=[d4w_xyuu; d4w_yyuu; d4w_yzuu; d4w_yuuu; i*d4w_yuuu; 0];
dCVar4(:,12)=[d4w_xzuu; d4w_yzuu; d4w_zzuu; d4w_zuuu; i*d4w_zuuu; 0];
dCVar4(:,13)=[d4w_xuuu; d4w_yuuu; d4w_zuuu; d4w_uuuu; i*d4w_uuuu; 0];
MVar4 = [u, w, dw_x, dw_y, dw_z, dw_u,...
    d2w_xx, d2w_xy, d2w_xz, d2w_yy, d2w_yz, d2w_zz,...
    d2w_xu, d2w_yu, d2w_zu, d2w_uu,...
    d3w_xxu, d3w_xyu, d3w_xzu, d3w_yyu, d3w_yzu, d3w_zzu,...
    d3w_xuu, d3w_yuu, d3w_zuu, d3w_uuu,...
    d4w_xxuu, d4w_xyuu, d4w_xzuu, d4w_yyuu, d4w_yzuu, d4w_zzuu,...
    d4w_xuuu, d4w_yuuu, d4w_zuuu, d4w_uuuu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableXYZ = {'x', 'y', 'z', 'u1', 'u2', 'gamma'};
dRiem_Gamma.x = sym('dRiem_Gamma_x',[6,6,6]);
dRiem_Gamma.y = sym('dRiem_Gamma_y',[6,6,6]);
dRiem_Gamma.z = sym('dRiem_Gamma_z',[6,6,6]);
dRiem_Gamma.u1 = sym('dRiem_Gamma_u1',[6,6,6]);
dRiem_Gamma.u2 = sym('dRiem_Gamma_u2',[6,6,6]);
dRiem_Gamma.gamma = sym('dRiem_Gamma_gamma',[6,6,6]);
% dRiem_Gamma.(variableXYZ{j})(m,n,k) = dRiem_Gamma(m,n,k)/dx_j.
ListRiemGamma = [];
for m=1:6
    for n=1:6
        for k=1:6
            ListRiemGamma = [ListRiemGamma;m,n,k];
        end
    end
end
% size(ListRiemGamma) = [216,3];
for count = 1:216
    m1 = ListRiemGamma(count,1);
    n1 = ListRiemGamma(count,2);
    k1 = ListRiemGamma(count,3);
    tempFun = Riem_Gamma(m1,n1,k1);
    tempDiffFun = dGamma_Feff_twistor_CR_funW(tempFun,RVar4,CVar4,dRVar4,dCVar4);
    for j=1:6
        dRiem_Gamma.(variableXYZ{j})(m1,n1,k1)=tempDiffFun(j);
    end
    clear tempFun
    clear tempDiffFun
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear complex
save('New_twistor_CR_funW_differentiate_Riem_Gamma_data.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    