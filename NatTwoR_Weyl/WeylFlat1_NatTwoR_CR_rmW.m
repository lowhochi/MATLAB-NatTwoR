load('DataWeyl2_NatTwoR_CR_rmW.mat');
% continue the previous Weyl6: M is flat.
clearvars wSet symSetWeylTwo %obsolete sets
% use wSetTwo variableSetW2
variable_NatTwoR_WeylFlat1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSetW2 = symvar(WeylTwo); %length = 69
subSetW2 = variableSetW2;
for j=1:69
    myVar = variableSetW2(j);
    myChar = char(myVar);
    if strfind(myChar,'dw')==1
       subSetW2(j)=myVar;
   elseif strfind(myChar,'d2w')==1
       subSetW2(j)=myVar;
   elseif strfind(myChar,'d3w')==1
       subSetW2(j)=myVar;
   elseif strfind(myChar,'d4w')==1
       subSetW2(j)=myVar;
   elseif myVar==u
       subSetW2(j)=myVar;
   elseif myVar==w
       subSetW2(j)=myVar;
   else
       subSetW2(j)=0;
       continue
    end
end
WeylFlat = sym('WeylFlat',[120,5]);
for j=1:120
    WeylFlat(j,1) = WeylTwo(j,1);
    WeylFlat(j,2) = WeylTwo(j,2);
    WeylFlat(j,3) = WeylTwo(j,3);
    WeylFlat(j,4)= WeylTwo(j,4);
    temp = WeylTwo(j,5);
    temp = subs(temp, variableSetW2, subSetW2);
    WeylFlat(j,5) = temp;
end
clearvars j temp myVar myChar
% checkArrayEqual(symvar(WeylFlat),wSetTwo)==1
Wflat1214 = WeylFlat(3,5);
Wflat1235 = WeylFlat(11,5);
Wflat1535 = WeylFlat(50,5);
Wflat1545 = WeylFlat(52,5);
Wflat1212 = WeylFlat(1,5);
Wflat1215 = WeylFlat(4,5);
Wflat1515 = WeylFlat(43,5);
Wflat1525 = WeylFlat(47,5);
WflatSet = [Wflat1214, Wflat1235, Wflat1535, Wflat1545,...
    Wflat1212, Wflat1215, Wflat1515, Wflat1525];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSet1 = [w, dw_u, dw_Mu, dw_conjMu, dw_vnormv, d2w_uu];
variableSet2 = [d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_MuMu,...
    d2w_MuconjMu, d2w_Muvnormv, d2w_conjMuconjMu,...
    d2w_conjMuvnormv, d2w_vnormvvnormv];
variableSet3 = [d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv, d3w_uuu];
variableSet4 = [d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv,...
    d3w_uconjMuconjMu, d3w_uconjMuvnormv, d4w_uuuu];
variableSet5 = [d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuconjMuconjMu,...
    d4w_uuuMu, d4w_uuuconjMu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
wSub = lambda0 +lambda1*u+ K*u^2 -conj(lambda1)*u^3 +conj(lambda0)*u^4;
dwVec = df_main4_CR_funWv2(wSub ,CVarFlat1, derivativeDict, gamma);
dw_uSub = dwVec(4);
dw_MuSub = dwVec(1);
dw_conjMuSub = dwVec(2);
dw_vnormvSub = dwVec(3);
d2w_uVec=df_main4_CR_funWv2(dw_uSub ,CVarFlat1, derivativeDict, gamma);
d2w_uuSub = d2w_uVec(4);
d2w_uMuSub = d2w_uVec(1);
d2w_uconjMuSub = d2w_uVec(2);
d2w_uvnormvSub = d2w_uVec(3);
d2w_MuVec=df_main4_CR_funWv2(dw_MuSub ,CVarFlat1, derivativeDict, gamma); 
d2w_conjMuVec=df_main4_CR_funWv2(dw_conjMuSub ,CVarFlat1, derivativeDict, gamma);
d2w_vnormvVec=df_main4_CR_funWv2(dw_vnormvSub ,CVarFlat1, derivativeDict, gamma);
d2w_MuMuSub = d2w_MuVec(1);
d2w_MuconjMuSub = d2w_MuVec(2);
d2w_MuvnormvSub = d2w_MuVec(3);
d2w_conjMuconjMuSub = d2w_conjMuVec(2);
d2w_conjMuvnormvSub = d2w_conjMuVec(3);
d2w_vnormvvnormvSub = d2w_vnormvVec(3);

d3w_uuVec=df_main4_CR_funWv2(d2w_uuSub,CVarFlat1,derivativeDict,gamma);
d3w_uuuSub = d3w_uuVec(4);
d3w_uuMuSub = d3w_uuVec(1);
d3w_uuconjMuSub = d3w_uuVec(2);
d3w_uuvnormvSub = d3w_uuVec(3);
d3w_uMuVec=df_main4_CR_funWv2(d2w_uMuSub,CVarFlat1,derivativeDict,gamma);
d3w_uconjMuVec=df_main4_CR_funWv2(d2w_uconjMuSub,CVarFlat1,derivativeDict,gamma);
d3w_uMuMuSub = d3w_uMuVec(1);
d3w_uMuconjMuSub = d3w_uMuVec(2);
d3w_uMuvnormvSub = d3w_uMuVec(3);
d3w_uconjMuconjMuSub = d3w_uconjMuVec(2);
d3w_uconjMuvnormvSub = d3w_uconjMuVec(3);

d4w_uuuVec=df_main4_CR_funWv2(d3w_uuuSub,CVarFlat1,derivativeDict,gamma);
d4w_uuMuVec=df_main4_CR_funWv2(d3w_uuMuSub,CVarFlat1,derivativeDict,gamma);
d4w_uuconjMuVec=df_main4_CR_funWv2(d3w_uuconjMuSub,CVarFlat1,derivativeDict,gamma);
d4w_uuuuSub = d4w_uuuVec(4);
d4w_uuuMuSub = d4w_uuuVec(1);
d4w_uuuconjMuSub = d4w_uuuVec(2);
d4w_uuMuMuSub = d4w_uuMuVec(1);
d4w_uuMuconjMuSub = d4w_uuMuVec(2);
d4w_uuconjMuconjMuSub = d4w_uuconjMuVec(2);

for number=1:8
    temp = WflatSet(number);
    for k=1:10
        temp = subs(temp, conjKset(k), Kset(k));
    end
    for j=1:length(variableSet5)
        myVar = variableSet5(j);
        myChar = [char(myVar),'Sub'];
        eval(['temp=subs(temp,',char(myVar),',',myChar,');']);
    end
    for j=1:length(variableSet4)
        myVar = variableSet4(j);
        myChar = [char(myVar),'Sub'];
        eval(['temp=subs(temp,',char(myVar),',',myChar,');']);
    end
    for j=1:length(variableSet3)
        myVar = variableSet3(j);
        myChar = [char(myVar),'Sub'];
        eval(['temp=subs(temp,',char(myVar),',',myChar,');']);
    end
    for j=1:length(variableSet2)
        myVar = variableSet2(j);
        myChar = [char(myVar),'Sub'];
        eval(['temp=subs(temp,',char(myVar),',',myChar,');']);
    end
    for j=1:length(variableSet1)
        myVar = variableSet1(j);
        myChar = [char(myVar),'Sub'];
        eval(['temp=subs(temp,',char(myVar),',',myChar,');']);
    end
    WflatSet(number) = complex_simple3(temp, MVarFlat1);
end

clearvars number k j temp myVar myChar
save('DataWeylFlat1_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charSet = ["Wflat1214", "Wflat1235", "Wflat1535", "Wflat1545",...
%     "Wflat1212", "Wflat1215", "Wflat1515", "Wflat1525"];
% latexFile = fopen('latex_WeylFlat1.txt','w');
% for number=1:8
%     temp= WflatSet(number);
%     latexTemp = latex(temp);
%     fprintf(latexFile,'%s : \\\\ \n', charSet(number));
%     fprintf(latexFile, '%s\n', ' ');
%     fprintf(latexFile,'$ %s $ \\\\[0.1in] \n', latexTemp);
%     fprintf(latexFile, '%s\n', ' ');
% end
% fclose(latexFile);
