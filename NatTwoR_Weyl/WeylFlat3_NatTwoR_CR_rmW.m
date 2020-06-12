% WeylFlat3_NatTwoR_CR_rmW.m
load('DataWeylFlat1_NatTwoR_CR_rmW.mat');
clearvars variableSet1 variableSet2 variableSet3 variableSet4 variableSet5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WflatSet02 = WflatSet(5:8);
% WflatSet = [Wflat1214, Wflat1235, Wflat1535, Wflat1545,...
%     Wflat1212, Wflat1215, Wflat1515, Wflat1525];
syms p y % p=x+i*z
syms phi % phi is real-valued
syms dphi_p dphi_y
dphi_conjp = conj(dphi_p);
syms d2phi_pp d2phi_pconjp d2phi_yy d2phi_py
d2phi_conjpy = conj(d2phi_py); 
d2phi_conjpconjp = conj(d2phi_pp);
syms d3phi_ppp d3phi_ppconjp 
d3phi_pconjpconjp = conj(d3phi_ppconjp);
d3phi_conjpconjpconjp = conj(d3phi_ppp);
syms d3phi_ppy d3phi_pconjpy d3phi_pyy d3phi_yyy
d3phi_conjpyy = conj(d3phi_pyy);
d3phi_conjpconjpy = conj(d3phi_ppy);
syms d4phi_pppp d4phi_pppconjp d4phi_ppconjpconjp
d4phi_pconjpconjpconjp = conj(d4phi_pppconjp);
d4phi_conjpconjpconjpconjp = conj(d4phi_pppp);
syms d4phi_pppy d4phi_ppconjpy 
d4phi_pconjpconjpy = conj(d4phi_ppconjpy);
d4phi_conjpconjpconjpy = conj(d4phi_pppy);
syms d4phi_ppyy d4phi_pconjpyy
d4phi_conjpconjpyy = conj(d4phi_ppyy);
syms d4phi_pyyy d4phi_yyyy
d4phi_conjpyyy = conj(d4phi_pyyy);
realVariable = [y, phi, dphi_y, d2phi_yy, d2phi_pconjp, d3phi_pconjpy,...
    d3phi_yyy, d4phi_ppconjpconjp, d4phi_pconjpyy, d4phi_yyyy];

variableSet = [lambda1, lambda0, K,...
    dlambda1_x, dlambda1_y, dlambda1_z, dlambda0_x, dlambda0_y, dlambda0_z,...
    dK_x, dK_y, dK_z, d2lambda1_xx, d2lambda1_xy, d2lambda1_xz, d2lambda1_yy,...
    d2lambda1_yz, d2lambda1_zz,  d2lambda0_xx, d2lambda0_xy, d2lambda0_xz,... 
    d2lambda0_yy, d2lambda0_yz, d2lambda0_zz,...
    d2K_xx, d2K_xy, d2K_xz, d2K_yy, d2K_yz, d2K_zz];

MVarFlat2 =[u, p, y, phi, dphi_p dphi_y,...
    d2phi_pp, d2phi_pconjp, d2phi_yy, d2phi_py,...
    d3phi_ppp, d3phi_ppconjp, d3phi_ppy, d3phi_pconjpy,...
    d3phi_pyy, d3phi_yyy, d4phi_pppp, d4phi_pppconjp,...
    d4phi_ppconjpconjp, d4phi_pppy, d4phi_ppconjpy,...
    d4phi_ppyy, d4phi_pconjpyy, d4phi_pyyy d4phi_yyyy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda1Sub = d2phi_py;
lambda0Sub = -1/2*d2phi_pp;
KSub = d2phi_pconjp-1/2*d2phi_yy;
% % % % %
dlambda1_xSub = d3phi_ppy + d3phi_pconjpy;
dlambda1_ySub = d3phi_pyy;
dlambda1_zSub = i*d3phi_ppy -i*d3phi_pconjpy;
dlambda0_xSub = -1/2*d3phi_ppp -1/2*d3phi_ppconjp;
dlambda0_ySub = -1/2*d3phi_ppy;
dlambda0_zSub = -i/2*d3phi_ppp +i/2*d3phi_ppconjp;
dK_xSub = d3phi_ppconjp + d3phi_pconjpconjp...
    -1/2*d3phi_pyy -1/2*d3phi_conjpyy;
dK_ySub = d3phi_pconjpy-1/2*d3phi_yyy;
dK_zSub = i*d3phi_ppconjp -i*d3phi_pconjpconjp...
    -i/2*d3phi_pyy +i/2*d3phi_conjpyy;
% % % % %
d2lambda1_xxSub = 2*d4phi_ppconjpy +d4phi_pconjpconjpy...
    + d4phi_pppy;
d2lambda1_xySub = d4phi_ppyy +d4phi_pconjpyy;
d2lambda1_xzSub = i*d4phi_pppy -i*d4phi_pconjpconjpy;
d2lambda1_yySub = d4phi_pyyy;
d2lambda1_yzSub =i*d4phi_ppyy -i*d4phi_pconjpyy;
d2lambda1_zzSub = 2*d4phi_ppconjpy -d4phi_pppy -d4phi_pconjpconjpy;
% % % % %
d2lambda0_xxSub = -d4phi_pppconjp -1/2*d4phi_pppp -1/2*d4phi_ppconjpconjp;
d2lambda0_xySub = -1/2*d4phi_pppy -1/2*d4phi_ppconjpy;
d2lambda0_xzSub = -i/2*d4phi_pppp +i/2*d4phi_ppconjpconjp;
d2lambda0_yySub = -1/2*d4phi_ppyy;
d2lambda0_yzSub = -i/2*d4phi_pppy +i/2*d4phi_ppconjpy;
d2lambda0_zzSub = 1/2*d4phi_pppp -d4phi_pppconjp +1/2*d4phi_ppconjpconjp;
% % % % %
d2K_xxSub = d4phi_pppconjp +2*d4phi_ppconjpconjp +d4phi_pconjpconjpconjp...
    -1/2*d4phi_ppyy -1/2*d4phi_conjpconjpyy -d4phi_pconjpyy;
d2K_xySub = d4phi_ppconjpy +d4phi_pconjpconjpy -1/2*d4phi_pyyy...
    -1/2*d4phi_conjpyyy;
d2K_xzSub = i*d4phi_pppconjp -i*d4phi_pconjpconjpconjp...
    -i/2*d4phi_ppyy +i/2*d4phi_conjpconjpyy;
d2K_yySub = d4phi_pconjpyy -1/2*d4phi_yyyy;
d2K_yzSub = i*d4phi_ppconjpy -i*d4phi_pconjpconjpy...
    -i/2*d4phi_pyyy +i/2*d4phi_conjpyyy;
d2K_zzSub = -d4phi_pppconjp +2*d4phi_ppconjpconjp...
    -d4phi_pconjpconjpconjp +1/2*d4phi_ppyy +1/2*d4phi_conjpconjpyy...
    -d4phi_pconjpyy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for number=1:4
    temp = WflatSet02(number);
    for j=1:length(variableSet)
        myVar = variableSet(j);
        myChar = [char(myVar),'Sub'];
        eval(['temp=subs(temp,',char(myVar),',',myChar,');']);
    end
    for k=1:10
        myVarTwo = realVariable(k);
        temp = subs(temp, conj(myVarTwo), myVarTwo);
    end
    WflatSet02(number)=complex_simple3(temp, MVarFlat2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars number temp j k myVar myChar myVarTwo
save('DataWeylFlat3_NatTwoR_CR_rmW.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
charSet = ["Wflat1212", "Wflat1215", "Wflat1515", "Wflat1525"];
latexFile = fopen('latex_WeylFlat3.txt','w');
for number=1:4
    temp = WflatSet02(number);
    fprintf(latexFile,'%s : \\\\ \n', charSet(number));
    fprintf(latexFile, '%s\n', ' ');
    [termVec, uVec] = coeffs(temp,[u, conj(u)]);
    for j=1:length(uVec)
        termTex = latex(termVec(j));
        uTex = latex(uVec(j));
        fprintf(latexFile,'term: \\spa $ %s $: \\\\ \n',uTex);
        fprintf(latexFile, '%s\n', ' ');
        fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n', termTex);
        fprintf(latexFile, '%s\n', ' ');
    end
end
fclose(latexFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

