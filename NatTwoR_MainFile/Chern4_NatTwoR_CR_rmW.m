% Chern4_NatTwoR_CR_rmW.m
load('DataChern2_NatTwoR_CR_rmW.mat');

variable_NatTwoR_Chern4_CR_rmW
% Apply exact model to compute ChernTwo(m,n,k,ll).
symSetChernTwo = symvar(ChernTwo);
% variableSetG =[];
% variableSetdG =[];
% variableSetd2G = [];
% wSet = [];
% for j=1:length(symSetChernTwo)
%     myVar = symSetChernTwo(j);
%     myString = string(myVar);
%     if strfind(myString,'G')==1
%         variableSetG = [variableSetG, myVar];
%     elseif strfind(myString,'dG')==1
%         variableSetdG = [variableSetdG, myVar];
%     elseif strfind(myString,'d2G')==1
%         variableSetd2G = [variableSetd2G, myVar];
%     else
%         wSet= [wSet,myVar];
%     end
% end

variableSetG = [G11_2, G11_3, G12_3, G22_1, G22_3, G23_1, G31_2, G33_1, G33_2];
variableSetdG = [dG11_2_by1, dG11_2_by2, dG11_3_by1, dG11_3_by2...
    dG11_3_by3, dG12_3_by1, dG12_3_by2, dG12_3_by3, dG22_1_by1,...
    dG22_1_by2, dG22_1_by3, dG22_3_by2, dG22_3_by3, dG23_1_by1,...
    dG23_1_by2, dG31_2_by1, dG31_2_by2, dG31_2_by3, dG33_1_by1,...
    dG33_1_by2, dG33_1_by3, dG33_2_by1, dG33_2_by2, dG33_2_by3];

variableSetd2G = [d2G11_2_by11, d2G11_2_by12, d2G11_2_by22, d2G11_3_by11,...
    d2G11_3_by12, d2G11_3_by13, d2G11_3_by22, d2G11_3_by23, d2G11_3_by33,...
    d2G12_3_by11, d2G12_3_by12, d2G12_3_by13, d2G12_3_by22, d2G12_3_by23,...
    d2G12_3_by33, d2G22_1_by11, d2G22_1_by12, d2G22_1_by13, d2G22_1_by22,...
    d2G22_1_by23, d2G22_1_by33, d2G22_3_by22, d2G22_3_by23, d2G22_3_by33,...
    d2G23_1_by11, d2G23_1_by12, d2G23_1_by13, d2G23_1_by22, d2G31_2_by11,...
    d2G31_2_by12, d2G31_2_by13, d2G31_2_by22, d2G31_2_by23, d2G31_2_by33,...
    d2G33_1_by11, d2G33_1_by12, d2G33_1_by13, d2G33_1_by22, d2G33_1_by23,...
    d2G33_1_by33, d2G33_2_by11, d2G33_2_by12, d2G33_2_by13, d2G33_2_by22,...
    d2G33_2_by23, d2G33_2_by33];

wSet = [d2w_conjMuconjMu, d2w_uconjMu, d2w_uu, dw_conjMu, dw_u, dw_vnormv, u, w];

% Torsion is zero when w = f.
% wSet = [d2w_conjMuconjMu, d2w_uconjMu, d2w_uu, dw_conjMu,...
%   dw_u, dw_vnormv, u, w]
f = -i/2*mu1*(mu3*G11_2 + mu1*G12_3 - mu2*G11_3)...
    -i/2*mu2*(-mu3*G22_1 + mu1*G22_3 + mu2*G23_1)...
    -i/2*mu3*(mu3*G31_2 - mu1*G33_2 + mu2*G33_1);
dfdu = complexdiff3(f,u,0);
d2fdu2 = complexdiff3(dfdu,u,0);
phiF = d2fdu2 - 6*conj(u)/(1+u*conj(u))*dfdu + 12*conj(u)^2/(1+u*conj(u))^2*f;

% Differentiation
dfVec = df_NatTwo_MuSet_CR_rmW(f, CVarG2, dCVarG2);
df_conjMu = dfVec(2);
df_vnormv = dfVec(3);
d2f_uVec = df_NatTwo_MuSet_CR_rmW(dfdu, CVarG2, dCVarG2);
d2f_uconjMu = d2f_uVec(2);
d2f_conjMuVec = df_NatTwo_MuSet_CR_rmW(df_conjMu, CVarG2, dCVarG2);
d2f_conjMuconjMu = d2f_conjMuVec(2);

% Simplification
df_conjMu = subs(df_conjMu, variableSetConjMu, subSetConjMu);
df_conjMu = subs(df_conjMu, variableSetBianchi1, subSetBianchi1);
df_vnormv = subs(df_vnormv, variableSetVnormv, subSetVnormv);
df_vnormv = subs(df_vnormv, variableSetBianchi1, subSetBianchi1);

d2f_uconjMu = subs(d2f_uconjMu, variableSetConjMu, subSetConjMu);
d2f_uconjMu = subs(d2f_uconjMu, variableSetBianchi1, subSetBianchi1);

d2f_conjMuconjMu = subs(d2f_conjMuconjMu, variableSetConjMuConjMu, subSetConjMuConjMu);
d2f_conjMuconjMu = subs(d2f_conjMuconjMu, variableSet0, subSet0);
d2f_conjMuconjMu = subs(d2f_conjMuconjMu, variableSetBianchi2, subSetBianchi2);
d2f_conjMuconjMu = subs(d2f_conjMuconjMu, variableSetBianchi1, subSetBianchi1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model: g0 = dx^2 + dy^2 + sin(y)^2*dz^2. 
% Only G33_2, dG33_2_by2, d2G33_2_by22 have non-zero values.
syms x1 x2 x3 real
xVec = [x1, x2, x3];
% metric g0 corresponding to dx, dy and dz.
g0 = [1, 0 ,0;
    0, 1, 0;
    0, 0, sin(x2)^2];
eVec1 = [1; 0; 0]; % e1 = 1*d/dx + 0*d/dy + 0*d/dz;
eVec2 = [0; 1; 0]; % e2 = 0*d/dx + 1*d/dy + 0*d/dz;
eVec3 = [0; 0; 1/sin(x2)]; % e3 = 0*d/dx + 0*d/dy + 1/sin(x2)*d/dz;
eMatrix = [eVec1, eVec2, eVec3];

% Find Gij_k and its derivatives automatically
lie_e1_e2 = sym('lie12',[3,1]);
lie_e2_e3 = sym('lie23',[3,1]);
lie_e1_e3 = sym('lie13',[3,1]);

for jj=1:3
    temp12 = 0;
    temp23 = 0;
    temp13 = 0;
    for ii=1:3
        temp12 = temp12 + eVec1(ii)*diff(eVec2(jj),xVec(ii))...
            - eVec2(ii)*diff(eVec1(jj),xVec(ii));
        temp23 = temp23 + eVec2(ii)*diff(eVec3(jj),xVec(ii))...
            - eVec3(ii)*diff(eVec2(jj),xVec(ii));
        temp13 = temp13 + eVec1(ii)*diff(eVec3(jj),xVec(ii))...
            - eVec3(ii)*diff(eVec1(jj),xVec(ii));
    end
    lie_e1_e2(jj) = temp12;
    lie_e2_e3(jj) = temp23;
    lie_e1_e3(jj) = temp13;
end 
% Rewrite lie_e1_e2 in the basis {e1,e2,e3}.
lie_e1_e2 = inv(eMatrix)*lie_e1_e2;
lie_e1_e3 = inv(eMatrix)*lie_e1_e3;
lie_e2_e3 = inv(eMatrix)*lie_e2_e3;

lieModel = cell(3,3);
% lieModel{i,j} = [e_i, e_j] on M.
lieModel{1,1} = zeros(3,1);
lieModel{2,2} = zeros(3,1);
lieModel{3,3} = zeros(3,1);

lieModel{1,2} = lie_e1_e2;
lieModel{1,3} = lie_e1_e3;
lieModel{2,3} = lie_e2_e3;
lieModel{2,1} = -lie_e1_e2;
lieModel{3,1} = -lie_e1_e3;
lieModel{3,2} = -lie_e2_e3;

Gs = sym('Gs',[3 3 3]);
% Gs(m,n,k) = Gs_{mn}^k with
% \nabla_{e_m}e_n = Gs_{mn}^k e_k. 
for m=1:3
    for n=1:3
        for k=1:3
            temp = -lieModel{n,k}(m) +lieModel{m,n}(k) +lieModel{k,m}(n);
            Gs(m,n,k)= 1/2*temp;     
        end
    end
end
%variableSetG = [G11_2, G11_3, G12_3, G22_1, G22_3, G23_1, G31_2, G33_1, G33_2];
subSetG = [Gs(1,1,2), Gs(1,1,3), Gs(1,2,3), Gs(2,2,1), Gs(2,2,3), Gs(2,3,1),...
    Gs(3,1,2), Gs(3,3,1), Gs(3,3,2)];

dGs = sym('dGs',[3 3 3 3]);
% dGs(m,n,k,p) = D(Gs_{mn}^k)) by e_p.
for m=1:3
    for n=1:3
        for k=1:3
            myTerm = Gs(m,n,k);
            myDeriv = [diff(myTerm,x1),diff(myTerm,x2),diff(myTerm,x3)];
            dGs(m,n,k,1) = myDeriv*eVec1;
            dGs(m,n,k,2) = myDeriv*eVec2;
            dGs(m,n,k,3) = myDeriv*eVec3;
        end
    end
end

% variableSetdG = [dG11_2_by1, dG11_2_by2, dG11_3_by1, dG11_3_by2...
%     dG11_3_by3, dG12_3_by1, dG12_3_by2, dG12_3_by3, dG22_1_by1,...
%     dG22_1_by2, dG22_1_by3, dG22_3_by2, dG22_3_by3, dG23_1_by1,...
%     dG23_1_by2, dG31_2_by1, dG31_2_by2, dG31_2_by3, dG33_1_by1,...
%     dG33_1_by2, dG33_1_by3, dG33_2_by1, dG33_2_by2, dG33_2_by3];
subSetdG = [dGs(1,1,2,1), dGs(1,1,2,2), dGs(1,1,3,1), dGs(1,1,3,2),...
    dGs(1,1,3,3), dGs(1,2,3,1), dGs(1,2,3,2), dGs(1,2,3,3), dGs(2,2,1,1),...
    dGs(2,2,1,2), dGs(2,2,1,3), dGs(2,2,3,2), dGs(2,2,3,3), dGs(2,3,1,1),...
    dGs(2,3,1,2), dGs(3,1,2,1), dGs(3,1,2,2), dGs(3,1,2,3), dGs(3,3,1,1),...
    dGs(3,3,1,2), dGs(3,3,1,3), dGs(3,3,2,1), dGs(3,3,2,2), dGs(3,3,2,3)];

d2Gs = sym('d2Gs',[3 3 3 3 3]);
% dGs(m,n,k,p,q) = D2(Gs_{mn}^k)) by e_p by e_q.
for m=1:3
    for n=1:3
        for k=1:3
            for p=1:3
                myTerm = dGs(m,n,k,p);
                myDeriv = [diff(myTerm,x1),diff(myTerm,x2),diff(myTerm,x3)];
                d2Gs(m,n,k,p,1) = myDeriv*eVec1;
                d2Gs(m,n,k,p,2) = myDeriv*eVec2;
                d2Gs(m,n,k,p,3) = myDeriv*eVec3;
            end
        end
    end
end

% variableSetd2G = [d2G11_2_by11, d2G11_2_by12, d2G11_2_by22, d2G11_3_by11,...
%     d2G11_3_by12, d2G11_3_by13, d2G11_3_by22, d2G11_3_by23, d2G11_3_by33,...
%     d2G12_3_by11, d2G12_3_by12, d2G12_3_by13, d2G12_3_by22, d2G12_3_by23,...
%     d2G12_3_by33, d2G22_1_by11, d2G22_1_by12, d2G22_1_by13, d2G22_1_by22,...
%     d2G22_1_by23, d2G22_1_by33, d2G22_3_by22, d2G22_3_by23, d2G22_3_by33,...
%     d2G23_1_by11, d2G23_1_by12, d2G23_1_by13, d2G23_1_by22, d2G31_2_by11,...
%     d2G31_2_by12, d2G31_2_by13, d2G31_2_by22, d2G31_2_by23, d2G31_2_by33,...
%     d2G33_1_by11, d2G33_1_by12, d2G33_1_by13, d2G33_1_by22, d2G33_1_by23,...
%     d2G33_1_by33, d2G33_2_by11, d2G33_2_by12, d2G33_2_by13, d2G33_2_by22,...
%     d2G33_2_by23, d2G33_2_by33];
subSetd2G = sym('subSetd2G',[1,length(variableSetd2G)]);
for j=1:length(variableSetd2G)
    myChar = char(variableSetd2G(j));
    m = str2double(myChar(4));
    n = str2double(myChar(5));
    k = str2double(myChar(7));
    p = str2double(myChar(11));
    q = str2double(myChar(12));
    subSetd2G(j) = d2Gs(m,n,k,p,q);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put w = f(xi,u) in the following context.
wSet = [d2w_conjMuconjMu, d2w_uconjMu, d2w_uu, dw_conjMu,...
  dw_u, dw_vnormv, u, w];
wSub = subs(f, variableSetG, subSetG);
dw_uSub = subs(dfdu, variableSetG, subSetG);
d2w_uuSub = subs(d2fdu2, variableSetG, subSetG);

dw_conjMuSub = subs(df_conjMu, variableSetdG, subSetdG);
dw_conjMuSub = subs(dw_conjMuSub, variableSetG, subSetG);
dw_vnormvSub = subs(df_vnormv, variableSetdG, subSetdG);
dw_vnormvSub = subs(dw_vnormvSub, variableSetG, subSetG);

d2w_uconjMuSub = subs(d2f_uconjMu, variableSetdG, subSetdG);
d2w_uconjMuSub = subs(d2w_uconjMuSub, variableSetG, subSetG);

d2w_conjMuconjMuSub = subs(d2f_conjMuconjMu, variableSetd2G, subSetd2G); 
d2w_conjMuconjMuSub = subs(d2w_conjMuconjMuSub, variableSetdG, subSetdG);
d2w_conjMuconjMuSub =  subs(d2w_conjMuconjMuSub, variableSetG, subSetG);

wSub = complex_simple3(wSub,u);
dw_uSub = complex_simple3(dw_uSub,u);
d2w_uuSub = complex_simple3(d2w_uuSub,u);
dw_conjMuSub = complex_simple3(dw_conjMuSub,u);
d2w_uconjMuSub = complex_simple3(d2w_uconjMuSub,u);
d2w_conjMuconjMuSub = complex_simple3(d2w_conjMuconjMuSub,u);

subSetW = [d2w_conjMuconjMuSub, d2w_uconjMuSub, d2w_uuSub, dw_conjMuSub,...
   dw_uSub, dw_vnormvSub, u, wSub];

% Substitution to ChernTwo tensor.
ChernModel = sym('ChernModel',[2 2 2 2]);

for j=1:16
    m = indexChern(j,1);
    n = indexChern(j,2);
    k = indexChern(j,3);
    ll = indexChern(j,4);
    temp = ChernTwo(m,n,k,ll);
    temp = subs(temp, wSet, subSetW);
    temp = subs(temp, variableSetd2G, subSetd2G);
    temp = subs(temp, variableSetdG, subSetdG);
    temp = subs(temp, variableSetG, subSetG);
    temp = complex_simple3(temp,[u]);
    ChernModel(m,n,k,ll) = temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars j myVar myChar m n k ll temp p q myTerm myDeriv ii jj





