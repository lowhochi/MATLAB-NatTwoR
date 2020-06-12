% flatM1_rewrite_NatTwoR_CR_rmW.m
% flatM1 and flatM2 refers to Chapter 3 content.

load('DataMain4_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RiemCurv0(i,j,k,l) = R_{ijk}^l
% RiemCurv(m,n,k,ll) = R_{mnkll}
% RicCurv(m,n) = Ric(u_m,u_n)
% "variableSet" = symSetRCurv
y = 1+u*conj(u);
subSetFlat = sym('subSet',[135,1]);
MVarFlat1 = [];
for j=1:135
   myVar = symSetRCurv(j);
   myChar = char(myVar);
   if strfind(myChar,'dw')==1
       subSetFlat(j)=myVar;
   elseif strfind(myChar,'d2w')==1
       subSetFlat(j)=myVar;
   elseif strfind(myChar,'d3w')==1
       subSetFlat(j)=myVar;
   elseif strfind(myChar,'drho')==1
       subSetFlat(j)=myVar;
   elseif strfind(myChar,'d2rho')==1
       subSetFlat(j)=myVar;
   elseif myVar==u
       subSetFlat(j)=myVar;
   elseif myVar==w
       subSetFlat(j)=myVar;
   elseif myVar==rho
       subSetFlat(j)=myVar; 
   else
       subSetFlat(j)=0;
       continue
   end
   MVarFlat1 = union(MVarFlat1,myVar);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(1): Find RicCurv when M is flat.
% Find RicSharp(m,n) = R_m^n = Ric#(u_m) on u_n component.
F0_Flat = sym('F0_Flat',[6 6]);
F0Inv_Flat = sym('F0Inv_Flat',[6 6]);
RicCurvFlat = sym('RicCurvFlat',[6 6]);
for m=1:6
    for n=1:6
        F0_Flat(m,n) = subs(F0(m,n), symSetRCurv, subSetFlat);
        F0Inv_Flat(m,n) = subs(F0Inv(m,n), symSetRCurv, subSetFlat);
        RicCurvFlat(m,n) = subs(RicCurv(m,n), symSetRCurv, subSetFlat);
    end
end

RicSharp = sym('RicSharpFlat',[6 6]);
for m=1:6
    for n=1:6
        temp = 0;
        for k=1:6
        temp = temp + RicCurvFlat(m,k)*F0Inv_Flat(k,n);
        end
        RicSharp(m,n) = complex_simple3(temp, MVarFlat1);
    end
end
clearvars m n temp j myVar myChar

variable_flatM1_CR_rmW

d3w_uuuZ = -i*(drho_u + 6*i*conj(u)/y*d2w_uu - 18*i*conj(u)^2/y^2*dw_u...
    +24*i*conj(u)^3/y^3*w -6*i/y^2*conj(dw_u) +24*i*u/y^3*conj(w));
d3w_uuconjMuZ = conj(d2GH11_1_conjuMu) + 2*conj(u)/y*conj(dGH11_1_Mu)...
    -2*conj(u)^2/y^2*conj(GH11_2);
d3w_uuMuZ = conj(d2GH11_1_conjuconjMu) +2*conj(u)/y*d2w_uMu...
    -2*conj(u)^2/y^2*dw_Mu;
d3w_uuvnormvZ = conj(d2GH11_1_conjuvnormv) +2*conj(u)/y*d2w_uvnormv...
    -2*conj(u)^2/y^2*dw_vnormv;
d2w_uMuZ = conj(dGH11_1_conjMu) + 2*conj(u)/y*dw_Mu;
d2w_uconjMuZ = conj(dGH11_1_Mu) + 2*conj(u)/y*dw_conjMu;
d2w_uvnormvZ = conj(dGH11_1_vnormv) + 2*conj(u)/y*dw_vnormv;
d2w_uuZ = phiW + 6*conj(u)/y*dw_u - 12*conj(u)^2/y^2*w;
dw_conjMuZ = -conj(GH11_2);

for m=1:6
    for n=1:6
    temp = RicSharp(m,n);

    temp = subs(temp, [d3w_uuu, d3w_uuconjMu, d3w_uuMu, d3w_uuvnormv],...
        [d3w_uuuZ, d3w_uuconjMuZ, d3w_uuMuZ, d3w_uuvnormvZ]);
    temp = subs(temp, [d2w_uMu, d2w_uconjMu, d2w_uvnormv], [d2w_uMuZ, d2w_uconjMuZ, d2w_uvnormvZ]);
    temp = subs(temp, [d2w_uu, dw_conjMu], [d2w_uuZ, dw_conjMuZ]);
    
    temp = subs(temp, conj(phiW), phiW+i*rho);
    temp = subs(temp, dw_u, conj(GH11_1)+2*conj(u)/y*w);
    temp = subs(temp, w*conj(w), -y^2*GT01_2);
    temp = subs(temp, u*conj(u), Y-1);
    temp = complex_simple3(temp, MVarFlat12);
    temp = subs(temp, u*conj(u), Y-1);
    temp = complex_simple3(temp, MVarFlat12);
    RicSharp(m,n) = temp;
    end
end
clear m n temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(2): Find RGamma when M is flat
RGammaFlat = sym('RGammaFlat',[6 6 6]);

for m=1:6
    for n=1:6
        for k=1:6
            temp = RGamma(m,n,k);
            temp = subs(temp, symSetRCurv, subSetFlat);
            temp = subs(temp, [d3w_uuu, d3w_uuconjMu, d3w_uuMu, d3w_uuvnormv],...
                [d3w_uuuZ, d3w_uuconjMuZ, d3w_uuMuZ, d3w_uuvnormvZ]);
            temp = subs(temp, [d2w_uMu, d2w_uconjMu, d2w_uvnormv],...
                [d2w_uMuZ, d2w_uconjMuZ, d2w_uvnormvZ]);
            temp = subs(temp, [d2w_uu, dw_conjMu], [d2w_uuZ, dw_conjMuZ]);
    
            temp = subs(temp, conj(phiW), phiW+i*rho);
            temp = subs(temp, dw_u, conj(GH11_1)+2*conj(u)/y*w);
            temp = subs(temp, w*conj(w), -y^2*GT01_2);
%             temp = subs(temp, u*conj(u), Y-1);
%             temp = complex_simple3(temp, MVarFlat12);
%             temp = subs(temp, u*conj(u), Y-1);
            temp = complex_simple3(temp, MVarFlat12);
            
            RGammaFlat(m,n,k) = temp;
        end
    end
end
clear m n k temp
save('DataFlatM1_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%