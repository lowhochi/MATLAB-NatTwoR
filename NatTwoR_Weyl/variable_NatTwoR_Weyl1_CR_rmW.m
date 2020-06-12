% variable_NatTwoR_Weyl1_CR_rmW.m

% d3aV_uXX
syms d3aV_uconjuconju d3aV_uuconju d3aV_uuu
syms d3aV_uMuMu d3aV_uMuconjMu d3aV_uMuvnormv d3aV_uconjMuconjMu 
syms d3aV_uconjMuvnormv d3aV_uvnormvvnormv
syms d3aV_uuMu d3aV_uuconjMu d3aV_uuvnormv 
syms d3aV_uconjuMu d3aV_uconjuconjMu d3aV_uconjuvnormv
d2aV_uRow = [d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv];
d3aV_uconjMuMu = d3aV_uMuconjMu + d2aV_uRow*lieMu{1,2};
d3aV_uvnormvMu = d3aV_uMuvnormv + d2aV_uRow*lieMu{1,3};
d3aV_uvnormvconjMu = d3aV_uconjMuvnormv + d2aV_uRow*lieMu{2,3};
d3aV_uMuu = d3aV_uuMu + d2aV_uRow*lieMu{4,1};
d3aV_uconjMuu = d3aV_uuconjMu + d2aV_uRow*lieMu{4,2};
d3aV_uvnormvu = d3aV_uuvnormv + d2aV_uRow*lieMu{4,3};
d3aV_uMuconju = d3aV_uconjuMu + d2aV_uRow*lieMu{5,1};
d3aV_uconjMuconju = d3aV_uconjuconjMu + d2aV_uRow*lieMu{5,2};
d3aV_uvnormvconju = d3aV_uconjuvnormv + d2aV_uRow*lieMu{5,3}; 

% d3aV_conjuXX
syms d3aV_conjuMuMu d3aV_conjuMuconjMu d3aV_conjuMuvnormv
syms d3aV_conjuconjMuconjMu d3aV_conjuconjMuvnormv d3aV_conjuvnormvvnormv
%d3aV_uconjuMu d3aV_uconjuconjMu d3aV_uconjuvnormv
syms d3aV_conjuconjuMu d3aV_conjuconjuconjMu d3aV_conjuconjuvnormv
% d3aV_uuconju d3aV_uconjuconju
syms d3aV_conjuconjuconju
d2aV_conjuRow = [d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv];

d3aV_conjuconjMuMu = d3aV_conjuMuconjMu + d2aV_conjuRow*lieMu{1,2};
d3aV_conjuvnormvMu = d3aV_conjuMuvnormv + d2aV_conjuRow*lieMu{1,3};
d3aV_conjuvnormvconjMu = d3aV_conjuconjMuvnormv + d2aV_conjuRow*lieMu{2,3};
d3aV_conjuMuu = d3aV_uconjuMu + d2aV_conjuRow*lieMu{4,1};
d3aV_conjuconjMuu = d3aV_uconjuconjMu + d2aV_conjuRow*lieMu{4,2};
d3aV_conjuvnormvu = d3aV_uconjuvnormv + d2aV_conjuRow*lieMu{4,3};
d3aV_conjuMuconju = d3aV_conjuconjuMu + d2aV_conjuRow*lieMu{5,1};
d3aV_conjuconjMuconju = d3aV_conjuconjuconjMu + d2aV_conjuRow*lieMu{5,2};
d3aV_conjuvnormvconju = d3aV_conjuconjuvnormv + d2aV_conjuRow*lieMu{5,3}; 

% d4aV_uconjuXX
syms d4aV_uconjuMuMu d4aV_uconjuMuconjMu d4aV_uconjuMuvnormv
syms d4aV_uconjuconjMuconjMu d4aV_uconjuconjMuvnormv d4aV_uconjuvnormvvnormv
syms d4aV_uuconjuMu d4aV_uuconjuconjMu d4aV_uuconjuvnormv
syms d4aV_uconjuconjuMu d4aV_uconjuconjuconjMu d4aV_uconjuconjuvnormv
syms d4aV_uuuconju d4aV_uuconjuconju d4aV_uconjuconjuconju

d3aV_uconjuRow = [d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv];
d4aV_uconjuconjMuMu = d4aV_uconjuMuconjMu + d3aV_uconjuRow*lieMu{1,2};
d4aV_uconjuvnormvMu = d4aV_uconjuMuvnormv + d3aV_uconjuRow*lieMu{1,3};
d4aV_uconjuvnormvconjMu = d4aV_uconjuconjMuvnormv + d3aV_uconjuRow*lieMu{2,3};
d4aV_uconjuMuu = d4aV_uuconjuMu + d3aV_uconjuRow*lieMu{4,1};
d4aV_uconjuconjMuu = d4aV_uuconjuconjMu + d3aV_uconjuRow*lieMu{4,2};
d4aV_uconjuvnormvu = d4aV_uuconjuvnormv + d3aV_uconjuRow*lieMu{4,3};
d4aV_uconjuMuconju = d4aV_uconjuconjuMu + d3aV_uconjuRow*lieMu{5,1};
d4aV_uconjuconjMuconju = d4aV_uconjuconjuconjMu + d3aV_uconjuRow*lieMu{5,2};
d4aV_uconjuvnormvconju = d4aV_uconjuconjuvnormv + d3aV_uconjuRow*lieMu{5,3};

% d4w_uuXX
syms d4w_uuMuMu d4w_uuMuconjMu d4w_uuMuvnormv d4w_uuconjMuconjMu
syms d4w_uuconjMuvnormv d4w_uuvnormvvnormv
syms d4w_uuuMu d4w_uuuconjMu d4w_uuuvnormv d4w_uuuu

d3w_uuRow =[d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv];
d4w_uuconjMuMu = d4w_uuMuconjMu + d3w_uuRow*lieMu{1,2};
d4w_uuvnormvMu = d4w_uuMuvnormv + d3w_uuRow*lieMu{1,3};
d4w_uuvnormvconjMu = d4w_uuconjMuvnormv + d3w_uuRow*lieMu{2,3};
d4w_uuMuu = d4w_uuuMu + d3w_uuRow*lieMu{4,1};
d4w_uuconjMuu = d4w_uuuconjMu + d3w_uuRow*lieMu{4,2};
d4w_uuvnormvu = d4w_uuuvnormv + d3w_uuRow*lieMu{4,3};
d4w_uuconjMuconju = d3w_uuRow*lieMu{5,2};
d4w_uuvnormvconju = d3w_uuRow*lieMu{5,3};

% dtheta
syms dtheta_Mu dtheta_vnormv
dtheta_conjMu = conj(dtheta_Mu);
% d2theta
syms d2theta_MuMu d2theta_MuconjMu d2theta_Muvnormv
syms d2theta_vnormvvnormv

d2theta_conjMuconjMu = conj(d2theta_MuMu); %theta is real
d2theta_conjMuMu = conj(d2theta_MuconjMu); %theta is real
d2theta_conjMuvnormv = conj(d2theta_Muvnormv); %theta is real

dthetaRow = [dtheta_Mu, dtheta_conjMu, dtheta_vnormv];
d2theta_vnormvMu = d2theta_Muvnormv + dthetaRow*lieMu{1,3}; %it has T4 aV bV h11
d2theta_vnormvconjMu = conj(d2theta_vnormvMu);  %theta is real
d2theta_Muu = dthetaRow*lieMu{4,1};
d2theta_conjMuconju = conj(d2theta_Muu); %theta is real
d2theta_vnormvu = dthetaRow*lieMu{4,3};
d2theta_vnormvconju = conj(d2theta_vnormvu);  %theta is real

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarWeyl = [u, aV, daV_u, daV_conju, d2aV_uu, d2aV_uconju,...
    w, dw_u, d2w_uu, theta,...
    daV_Mu, daV_conjMu, daV_vnormv, dtheta_Mu,...
    d2aV_uMu, d2aV_uconjMu, d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuconju,...
    d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv, d3aV_uuconju,...
    d3aV_uconjuconju,...
    dw_Mu, dw_conjMu, d2w_uMu, d2w_uconjMu, d3w_uuMu, d3w_uuconjMu, d3w_uuu];

myNumber = length(CVarWeyl);
dCVarWeyl = sym('dCVar',[5,myNumber]);

dCVarWeyl(:,1)=[0; 0; 0; 1; 0]; %u
dCVarWeyl(:,2)=[daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarWeyl(:,3)=[d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv; d2aV_uu; d2aV_uconju]; %daV_u
dCVarWeyl(:,4)=[d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;...
    d2aV_uconju; d2aV_conjuconju]; %daV_conju
dCVarWeyl(:,5)=[d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv; d3aV_uuu;...
    d3aV_uuconju]; %d2aV_uu
dCVarWeyl(:,6)=[d3aV_uconjuMu; d3aV_uconjuconjMu; d3aV_uconjuvnormv;...
    d3aV_uuconju; d3aV_uconjuconju];%d2aV_uconju
dCVarWeyl(:,7)=[dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0]; %w
dCVarWeyl(:,8)=[d2w_uMu; d2w_uconjMu; d2w_uvnormv; d2w_uu; 0]; %dw_u
dCVarWeyl(:,9)=[d3w_uuMu; d3w_uuconjMu; d3w_uuvnormv; d3w_uuu; 0]; %d2w_uu
dCVarWeyl(:,10)=[dtheta_Mu; dtheta_conjMu; dtheta_vnormv;
    0; 0]; %theta

dCVarWeyl(:,11)=[d2aV_MuMu; d2aV_MuconjMu; d2aV_Muvnormv; d2aV_Muu; 
    d2aV_Muconju]; %daV_Mu
dCVarWeyl(:,12)=[d2aV_conjMuMu; d2aV_conjMuconjMu; d2aV_conjMuvnormv;
    d2aV_conjMuu; d2aV_conjMuconju]; %daV_conjMu
dCVarWeyl(:,13)=[d2aV_vnormvMu; d2aV_vnormvconjMu; d2aV_vnormvvnormv; 
    d2aV_vnormvu; d2aV_vnormvconju]; %daV_vnormv
dCVarWeyl(:,14)=[d2theta_MuMu; d2theta_MuconjMu; d2theta_Muvnormv;
    d2theta_Muu; 0]; %dtheta_Mu

dCVarWeyl(:,15)=[d3aV_uMuMu; d3aV_uMuconjMu; d3aV_uMuvnormv; d3aV_uMuu;
    d3aV_uMuconju]; %d2aV_uMu
dCVarWeyl(:,16)=[d3aV_uconjMuMu; d3aV_uconjMuconjMu; d3aV_uconjMuvnormv; 
    d3aV_uconjMuu; d3aV_uconjMuconju]; %d2aV_uconjMu
dCVarWeyl(:,17)=[d3aV_conjuMuMu; d3aV_conjuMuconjMu; d3aV_conjuMuvnormv;
    d3aV_conjuMuu; d3aV_conjuMuconju]; %d2aV_conjuMu
dCVarWeyl(:,18)=[d3aV_conjuconjMuMu; d3aV_conjuconjMuconjMu; d3aV_conjuconjMuvnormv;
    d3aV_conjuconjMuu; d3aV_conjuconjMuconju]; %d2aV_conjuconjMu
dCVarWeyl(:,19)=[d3aV_conjuconjuMu; d3aV_conjuconjuconjMu; d3aV_conjuconjuvnormv;
    d3aV_uconjuconju; d3aV_conjuconjuconju]; %d2aV_conjuconju

dCVarWeyl(:,20)=[d4aV_uconjuMuMu; d4aV_uconjuMuconjMu; d4aV_uconjuMuvnormv;
    d4aV_uconjuMuu; d4aV_uconjuMuconju]; %d3aV_uconjuMu
dCVarWeyl(:,21)=[d4aV_uconjuconjMuMu; d4aV_uconjuconjMuconjMu;
    d4aV_uconjuconjMuvnormv; d4aV_uconjuconjMuu;
    d4aV_uconjuconjMuconju]; %d3aV_uconjuconjMu
dCVarWeyl(:,22)=[d4aV_uconjuvnormvMu; d4aV_uconjuvnormvconjMu;
    d4aV_uconjuvnormvvnormv; d4aV_uconjuvnormvu; d4aV_uconjuvnormvconju]; 
    %d3aV_uconjuvnormv
dCVarWeyl(:,23)=[d4aV_uuconjuMu; d4aV_uuconjuconjMu; d4aV_uuconjuvnormv;
    d4aV_uuuconju; d4aV_uuconjuconju];%d3aV_uuconju
dCVarWeyl(:,24)=[d4aV_uconjuconjuMu; d4aV_uconjuconjuconjMu;
    d4aV_uconjuconjuvnormv; d4aV_uuconjuconju; d4aV_uconjuconjuconju];
    %d3aV_uconjuconju
dCVarWeyl(:,25)=[d2w_MuMu; d2w_MuconjMu; d2w_Muvnormv; d2w_Muu; 0]; %dw_Mu
dCVarWeyl(:,26)=[d2w_conjMuMu; d2w_conjMuconjMu; d2w_conjMuvnormv; 
    d2w_conjMuu; d2w_conjMuconju]; %dw_conjMu
dCVarWeyl(:,27)=[d3w_uMuMu; d3w_uMuconjMu; d3w_uMuvnormv; d3w_uMuu;
    0]; %d2w_uMu
dCVarWeyl(:,28)=[d3w_uconjMuMu; d3w_uconjMuconjMu; d3w_uconjMuvnormv;
    d3w_uconjMuu; d3w_uconjMuconju]; %d2w_uconjMu

dCVarWeyl(:,29)=[d4w_uuMuMu; d4w_uuMuconjMu; d4w_uuMuvnormv;
    d4w_uuMuu; 0]; %d3w_uuMu
dCVarWeyl(:,30)=[d4w_uuconjMuMu; d4w_uuconjMuconjMu; 
    d4w_uuconjMuvnormv; d4w_uuconjMuu; d4w_uuconjMuconju]; %d3w_uuconjMu
dCVarWeyl(:,31)=[d4w_uuuMu; d4w_uuuconjMu; d4w_uuuvnormv; d4w_uuuu; 0]; %d3w_uuu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MVarWeyls1 = [wSet, aVSet];
MVarWeyls2 = [d3aV_uconjuconju, d3aV_uuconju, d3aV_uuu,...
    d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv, d3aV_uconjMuconjMu,...
    d3aV_uconjMuvnormv, d3aV_uvnormvvnormv, d3aV_uuMu, d3aV_uuconjMu,...
    d3aV_uuvnormv, d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv,...
    d3aV_conjuMuMu, d3aV_conjuMuconjMu, d3aV_conjuMuvnormv,...
    d3aV_conjuconjMuconjMu, d3aV_conjuconjMuvnormv, d3aV_conjuvnormvvnormv,...
    d3aV_conjuconjuMu, d3aV_conjuconjuconjMu, d3aV_conjuconjuvnormv,...
    d3aV_conjuconjuconju,...
    d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
    d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuvnormvvnormv,...
    d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuvnormv,...
    d4aV_uconjuconjuMu, d4aV_uconjuconjuconjMu, d4aV_uconjuconjuvnormv,...
    d4aV_uuuconju, d4aV_uuconjuconju, d4aV_uconjuconjuconju,...
    d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuMuvnormv, d4w_uuconjMuconjMu,...
    d4w_uuconjMuvnormv, d4w_uuvnormvvnormv,...
    d4w_uuuMu, d4w_uuuconjMu, d4w_uuuvnormv, d4w_uuuu,...
    dtheta_Mu, dtheta_conjMu, dtheta_vnormv,...
    d2theta_MuMu, d2theta_MuconjMu, d2theta_Muvnormv, d2theta_vnormvvnormv];

MVarWeylOne = union(MVarWeyls1, MVarWeyls2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countSet225 = [];
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                rowTemp = [m,n,k,ll];
                if (m<n)&&(k<ll)
                    countSet225 = [countSet225; rowTemp];
                end
            end
        end
    end
end
numberCount = zeros(225,2);
countIndex120 = []; % length(countSetTwo)=120
for j=1:225
    numberCount(j,1) = countSet225(j,1)*10 + countSet225(j,2);
    numberCount(j,2) = countSet225(j,3)*10 + countSet225(j,4);
    if numberCount(j,1)<=numberCount(j,2)
        countIndex120 = [countIndex120; countSet225(j,:)];
    end      
end
clearvars k numberCount countSet225 m n k ll rowTemp numberCount
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
