% example2_Part2_NatTwoR_CR_rmW.m

load('DataWeyl2_NatTwoR_CR_rmW.mat');
clearvars myChar myCharR
% Wf1212 = Weylf(1,5);
% Wf1215 = Weylf(4,5);
% Wf1515 = Weylf(43,5);
% Wf1525 = Weylf(47,5);
variable_eg2_NatTwoR_CR_rmW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSet1 = [aV, daV_u, daV_conju, daV_Mu, daV_conjMu, daV_vnormv];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% substitution on WeylTwo
WeyltTwo = sym('WeyltTwo',[120,5]);

for j=1:120
    for k=1:4
        WeyltTwo(j,k) = WeylTwo(j,k);
    end
    temp = WeylTwo(j,5); 
    % % % % %
    for m=1:length(variableSetf4)
        mychar = char(variableSetf4(m));
        mycharf = strrep(mychar,'w','f');
        eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
    end
    % % % % %
    for m=1:length(variableSetf3)
        mychar = char(variableSetf3(m));
        mycharf = strrep(mychar,'w','f');
        eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
    end
    % % % % %
    for m=1:length(variableSetf2)
        mychar = char(variableSetf2(m));
        mycharf = strrep(mychar,'w','f');
        eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
    end
    % % % % %
    for m=1:length(variableSetf1)
        mychar = char(variableSetf1(m));
        mycharf = strrep(mychar,'w','f');
        eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
    end
    % % % % %
    for m=1:length(variableSetT)
        mychar = char(variableSetT(m));
        mycharT = insertAfter(mychar,'a','s');
        eval(['temp=subs(temp,', mychar, ',', mycharT, ');']);
    end
    % % % % %
    for m=1:length(variableSet4)
        mychar = char(variableSet4(m));
        mychars = insertAfter(mychar,'V','s');
        eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
    end
    % % % % %
    for m=1:length(variableSet3)
        mychar = char(variableSet3(m));
        mychars = insertAfter(mychar,'V','s');
        eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
    end
    % % % % %
    for m=1:length(variableSet2)
        mychar = char(variableSet2(m));
        mychars = insertAfter(mychar,'V','s');
        eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
    end
    % % % % %
    for m=1:length(variableSet1)
        mychar = char(variableSet1(m));
        mychars = insertAfter(mychar,'V','s');
        eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
    end    
    WeyltTwo(j,5) = complex_simple3(temp, u);
end
clearvars m j temp k temp0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataTemp_Dec19.mat', 'WeyltTwo');

% for j=1:120
%     temp = WeyltTwo(j,5);
%     if temp~=0
%         disp(j);
%         disp(WeyltTwo(j,1:4));
%     end
% end
W1212 = WeyltTwo(1,5);
W1515 = WeyltTwo(43,5);

factor1212 = factor(W1212); % length = 7
