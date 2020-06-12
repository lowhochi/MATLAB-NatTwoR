% example_detd2f_being_zero.m
syms x y z 
syms f df_x df_y df_z
syms d2f_xx d2f_xy d2f_xz d2f_yy d2f_yz d2f_zz
syms d3f_xxx d3f_xxy d3f_xxz d3f_xyy d3f_xyz d3f_xzz
syms d3f_yyy d3f_yyz d3f_yzz d3f_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2f = [d2f_xx, d2f_xy, d2f_xz;
    d2f_xy, d2f_yy, d2f_yz;
    d2f_xz, d2f_yz, d2f_zz];
d3fRow = [d3f_xxx, d3f_xxy, d3f_xxz, d3f_xyy, d3f_xyz, d3f_xzz,...
    d3f_yyy, d3f_yyz, d3f_yzz, d3f_zzz];
MVarf = [x, y, z, f, df_x, df_y, df_z, ...
    d2f_xx, d2f_xy, d2f_xz, d2f_yy, d2f_yz, d2f_zz,...
    d3f_xxx, d3f_xxy, d3f_xxz, d3f_xyy, d3f_xyz, d3f_xzz, ...
    d3f_yyy, d3f_yyz, d3f_yzz, d3f_zzz];

CVar = [x, y, z, f, df_x, df_y, df_z,...
    d2f_xx, d2f_xy, d2f_xz, d2f_yy, d2f_yz, d2f_zz];

derivativeDict.x = [1; 0; 0];
derivativeDict.y = [0; 1; 0];
derivativeDict.z = [0; 0; 1];
derivativeDict.f = [df_x; df_y; df_z];
derivativeDict.df_x = [d2f_xx; d2f_xy; d2f_xz];
derivativeDict.df_y = [d2f_xy; d2f_yy; d2f_yz];
derivativeDict.df_z = [d2f_xz; d2f_yz; d2f_zz];
derivativeDict.d2f_xx = [d3f_xxx; d3f_xxy; d3f_xxz];
derivativeDict.d2f_xy = [d3f_xxy; d3f_xyy; d3f_xyz];
derivativeDict.d2f_xz = [d3f_xxz; d3f_xyz; d3f_xzz];
derivativeDict.d2f_yy = [d3f_xyy; d3f_yyy; d3f_yyz];
derivativeDict.d2f_yz = [d3f_xyz; d3f_yyz; d3f_yzz];
derivativeDict.d2f_zz = [d3f_xzz; d3f_yzz; d3f_zzz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume d2f_xx =/=0 and -d2f_xy^2+d2f_xx*d2f_yy =/= 0
step01 = myRowOperation(d2f, 2, d2f_xx, 1, -d2f_xy, MVarf);
step02 = myRowOperation(step01, 3, d2f_xx, 1, -d2f_xz, MVarf);
step03 = myRowOperation(step02, 3, d2f_xy^2 -d2f_xx*d2f_yy,...
    2, step02(3,2), MVarf);
step03(3,3) = 0; % % % % det(d2f)=0
step04 = myRowOperation(step03, 1, d2f_xy^2 -d2f_xx*d2f_yy,...
    2, d2f_xy, MVarf);

% aTerm = -(d2f_xy*d2f_yz - d2f_xz*d2f_yy)/(- d2f_xy^2 + d2f_xx*d2f_yy);
% bTerm = (d2f_xx*d2f_yz - d2f_xy*d2f_xz)/(- d2f_xy^2 + d2f_xx*d2f_yy);
aTerm = step04(1,3)/step04(1,1);
bTerm = step04(2,3)/step04(2,2);
aTerm = complex_simple3(aTerm, MVarf);
bTerm = complex_simple3(bTerm, MVarf);
% column: [daTerm_x, daTerm_y, daTerm_z, dbTerm_x, dbTerm_y, dbTerm_z];
% first derivative matrix of aTerm and bTerm
derivMatAB = [-d2f_xy, d2f_xx, 0, -d2f_yy, d2f_xy, 0; %2.1
    -d2f_xz, 0, d2f_xx, -d2f_yz, 0, d2f_xy; %2.2
    0, -d2f_xz, d2f_xy, 0, -d2f_yz, d2f_yy]; %2.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vvec = [-aTerm; -bTerm; 1];
d3f_byV = sym('d3f_byV',[3 3]);
for j=1:3
    for k=1:3
        d3fijVec = df_example_xyz(d2f(j,k),CVar,derivativeDict);
        temp = d3fijVec*Vvec;
        d3f_byV(j,k) = complex_simple3(temp, MVarf);
        clear temp
    end
end
prodVec01 = d3f_byV*[aTerm; bTerm; -1];
for k=1:3
    prodVec01(k) = expand(prodVec01(k));
    prodVec01(k) = complex_simple3(prodVec01(k), MVarf);
end
clearvars j k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aPterm = prodVec01(1);
syms C33 C13 C23
aPterm = subs(aPterm, d2f_xx*d2f_yy-d2f_xy^2, C33);
aPterm = complex_simple3(aPterm, [MVarf,C13,C23,C33]);
% [cfVec, fVec] = coeffs(aPterm, d3fRow); 

bPterm = prodVec01(2);
bPterm = subs(bPterm, d2f_xx*d2f_yy-d2f_xy^2, C33);
bPterm = complex_simple3(bPterm, [MVarf,C13,C23,C33]);
[cfVec02, fVec02] = coeffs(bPterm, d3fRow); 
for j=1:length(fVec02);
    disp(fVec02(j));
    temp = factor(cfVec02(j));
    disp(temp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





