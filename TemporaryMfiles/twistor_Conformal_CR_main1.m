% twistor_Conformal_CR_main1.m
% M is a flat space. alphaC = exp(2*f)*alpha
variable_Conformal_main1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in [e1, e2, e3, ddu, ddconju];
conjX1 = [u^2-1; 2*u; i*(u^2+1); w; 0];
X1 = [conj(u)^2-1; 2*conj(u); -i*(conj(u)^2+1); 0; conj(w)];
X2 = [0;0;0;1;0];
conjX2 = [0;0;0;0;1];
% difference in Tc and alphaC
df_conjX1 = df_Mu + w*df_u;
df_X1 = df_conjMu + conj(w)*df_conju;

Tc = exp(-2*f)*[v1normv; v2normv; v3normv; 0; 0]...
    +exp(-2*f)*(df_conju*X1 +df_u*conjX1 -df_conjX1*X2 - df_X1*conjX2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Webster metric gC in the basis (X1, conj_X1, X2, conj_X2, Tc)
hC = [0, -i*exp(2*f);
    i*exp(2*f), 0];
hCinv = [0, -i*exp(-2*f);
    i*exp(-2*f), 0];
gC = [0, 0, 0, -i*exp(2*f), 0;
    0, 0, i*exp(2*f), 0, 0;
    0, i*exp(2*f), 0, 0, 0;
    -i*exp(2*f), 0, 0, 0, 0;
    0, 0, 0, 0, 1];
gCinv = inv(gC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uvector = sym('Uvector',[5 5]); %in [Mu; conjMu; vnormv; u; conju];
Uvector(:,1) = [0; 1; 0; 0; conj(w)]; %X1
Uvector(:,2) = [1; 0; 0; w; 0]; %conjX1
Uvector(:,3) = [0; 0; 0; 1; 0]; %X2
Uvector(:,4) = [0; 0; 0; 0; 1]; %conjX2
Uvector(:,5) = exp(-2*f)*[0; 0; 1; 0; 0]...
    +exp(-2*f)*df_conju*[0; 1; 0; 0; conj(w)]...
    +exp(-2*f)*df_u*[1; 0; 0; w; 0]...    
    -exp(-2*f)*df_conjX1*[0; 0; 0; 1; 0]...
    -exp(-2*f)*df_X1*[0; 0; 0; 0; 1]; %Tc

d2f_X1Vec = df_Conformal_main1(df_X1,CVar1,derivativeDict);
d2f_X1X1 = d2f_X1Vec*Uvector(:,1);
d2f_X1conjX1 = d2f_X1Vec*Uvector(:,2);
d2f_conjX1conjX1 = conj(d2f_X1X1);
d2f_X1u = d2f_X1Vec(4);
d2f_X1conju = d2f_X1Vec(5);
d2f_conjX1u = conj(d2f_X1conju);

d2f_uX1 = [d2f_uMu,d2f_uconjMu,d2f_uvnormv,d2f_uu,d2f_uconju]*Uvector(:,1);
d2f_uconjX1 = [d2f_uMu,d2f_uconjMu,d2f_uvnormv,d2f_uu,d2f_uconju]*Uvector(:,2);
d2f_conjuconjX1 = conj(d2f_uX1);
% % % % %
lieOne = cell(5,5);
% Lie vectors are in (X1, conjX1, X2, conjX2, Tc)
% lieOne{1,:} = {0, [X1,conjX1], [X1,X2], [X1,conjX2], [X1,Tc]};
% lieOne{2,:} = {-[X1,conjX1], 0, [conjX1,X2], [conjX1,conjX2], [conjX1,Tc]};
% lieOne{3,:} = {-[X1,X2], -[conjX1,X2], 0, [X2,conjX2], [X2,Tc]};
% lieOne{4,:} = {-[X1,conjX2], -[conjX1,conjX2], -[X2,conjX2], 0, [conjX2,Tc]};
% lieOne{5,:} = {-[X1,Tc], -[conjX1,Tc], -[X2,Tc], -[conjX2,Tc], 0};
lie_X1_conjX1 = [0; 0; dw_conjMu; -conj(dw_conjMu); 0];
lie_X1_X2 = zeros(5,1);
lie_X1_conjX2 = [-2*u/(1+u*conj(u))+2*df_conju;
    2*df_u;
    -2*df_conjX1;
    -conj(dw_u)+2*u*conj(w)/(1+u*conj(u))-2*df_X1;
    -2*exp(2*f)];
lie_conjX1_X2 = [conj(lie_X1_conjX2(2)); 
    conj(lie_X1_conjX2(1));
    conj(lie_X1_conjX2(4));
    conj(lie_X1_conjX2(3));
    conj(lie_X1_conjX2(5))];
% % % % Involving Tc
lie_conjX1_Tc = -2*df_conjX1*[0; 0; 0; 0; 1]...
    +exp(-2*f)*[-w/(1+u*conj(u))^2; 0; -dw_vnormv; 
        w*conj(w)/(1+u*conj(u))^2; 0]...
    +exp(-2*f)*[d2f_conjuconjX1; d2f_uconjX1; 
        -d2f_conjX1conjX1; -d2f_X1conjX1; 0]...
    +exp(-2*f)*(-df_conju*lie_X1_conjX1-df_conjX1*lie_conjX1_X2);
% % % % %    
lie_X2_Tc = -2*df_u*[0; 0; 0; 0; 1]...
    +exp(-2*f)*[-1/(1+u*conj(u))^2; 0; 0; conj(w)/(1+u*conj(u))^2; 0]...
    +exp(-2*f)*[d2f_uconju; d2f_uu; -d2f_conjX1u; -d2f_X1u; 0]...
    -exp(-2*f)*df_u*lie_conjX1_X2;
% % % % %   
lie_X1_Tc = [conj(lie_conjX1_Tc(2));
    conj(lie_conjX1_Tc(1));
    conj(lie_conjX1_Tc(4));
    conj(lie_conjX1_Tc(3));
    conj(lie_conjX1_Tc(5))];
% % % % %   
lie_conjX2_Tc = [conj(lie_X2_Tc(2));
    conj(lie_X2_Tc(1));
    conj(lie_X2_Tc(4));
    conj(lie_X2_Tc(3));
    conj(lie_X2_Tc(5))];
% % % % % Row of X1
lieOne{1,1} = zeros(5,1);
lieOne{1,2} = lie_X1_conjX1;
lieOne{1,3} = zeros(5,1);
lieOne{1,4} = lie_X1_conjX2;
lieOne{1,5} = lie_X1_Tc;
% Row of conjX1
lieOne{2,1} = -lie_X1_conjX1;
lieOne{2,2} = zeros(5,1);
lieOne{2,3} = lie_conjX1_X2;
lieOne{2,4} = zeros(5,1);
lieOne{2,5} = lie_conjX1_Tc;
% Row of X2
lieOne{3,1} = zeros(5,1);
lieOne{3,2} = -lie_conjX1_X2;
lieOne{3,3} = zeros(5,1);
lieOne{3,4} = zeros(5,1);
lieOne{3,5} = lie_X2_Tc;
% Row of conjX2
lieOne{4,1} = -lie_X1_conjX2;
lieOne{4,2} = zeros(5,1);
lieOne{4,3} = zeros(5,1);
lieOne{4,4} = zeros(5,1);
lieOne{4,5} = lie_conjX2_Tc;
% Row of Tc
lieOne{5,1} = -lie_X1_Tc;
lieOne{5,2} = -lie_conjX1_Tc;
lieOne{5,3} = -lie_X2_Tc;
lieOne{5,4} = -lie_conjX2_Tc;
lieOne{5,5} = zeros(5,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dhC = cell(2,2);
% dhC{m,n} = in [Mu, conjMu, vnormv, u, conju]
dhC{1,1} = [0, 0, 0, 0, 0];
dhC{1,2} = -2*i*exp(2*f)*[df_Mu, df_conjMu, df_vnormv, df_u, df_conju];
dhC{2,1} = 2*i*exp(2*f)*[df_Mu, df_conjMu, df_vnormv, df_u, df_conju];
dhC{2,2} = [0, 0, 0, 0, 0];

GammaC.holo = sym('GammaC_holo_%d_%d_%d',[2 2 2]);
GammaC.antiholo = sym('GammaC_antiholo_%d_%d_%d',[2 2 2]);
GammaC.T = sym('GammaC_T_%d_%d', [2 2]);
for m=1:2
    for n=1:2
        for k=1:2
            temp1 = 0;
            temp2 = 0;
            for ll=1:2
                temp1 = temp1 + hCinv(ll,k)*(dhC{n,ll}*Uvector(:,2*m-1)...
                    -gC(2*n-1,:)*lieOne{2*m-1,2*ll});
                temp2 = temp2 + hCinv(ll,k)*(gC(2*ll,:)*lieOne{2*m, 2*n-1});
            end
            GammaC.holo(m,n,k) = temp1;
            GammaC.antiholo(m,n,k) = temp2;
            clearvars temp1 temp2
        end
    end
end
for n=1:2
    for k=1:2
        tempT = 0;
        for ll=1:2
            tempT = tempT + hCinv(ll,k)*(gC(2*ll,:)*lieOne{5,2*n-1});
        end
        GammaC.T(n,k) = tempT;
        clear tempT
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
realVariable = [f, df_vnormv, d2f_MuconjMu, d2f_vnormvvnormv, d2f_uconju];
conjRealVariable = sym('crV',[1,length(realVariable)]);
for j=1:length(realVariable)
    tempVar = realVariable(j);
    conjRealVariable(j) = conj(tempVar);
end

for n=1:2
    for k=1:2
        temp0 = GammaC.T(n,k);
        temp0 = subs(temp0, conjRealVariable, realVariable);
        GammaC.T(n,k) = complex_simple3(temp0, MVar1);
        for m=1:2
            temp1 = GammaC.holo(m,n,k);
            temp2 = GammaC.antiholo(m,n,k);
            temp1 = subs(temp1, conjRealVariable, realVariable);
            temp2 = subs(temp2, conjRealVariable, realVariable);
            GammaC.holo(m,n,k)= complex_simple3(temp1, MVar1);
            GammaC.antiholo(m,n,k)= complex_simple3(temp2, MVar1);
        end
    end
end
clearvars j n k ll m temp1 temp2 temp0 tempVar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Data_Conformal_CR_Main1.mat');