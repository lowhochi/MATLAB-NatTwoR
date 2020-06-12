function df = df_NatTwo_MuGamma_CR_rmW(f, CVar1, dCVar1, gamma)
% Output: df = [df/dMu, df/dconjMu, df/dvnormv, df/du, df/dconju, df/dgamma]
% CVar may contain both real and complex variables.
% dCVar1(1,:) = differentiate CVar1(:) by Mu;
% dCVar1(2,:) = differentiate CVar1(:) by conjMu;
% dCVar1(3,:) = differentiate CVar1(:) by vnormv;
% dCVar1(4,:) = differentiate CVar1(:) by u;
% dCVar1(5,:) = differentiate CVar1(:) by conju;
% dCVar1(6,:) = differentiate CVar1(:) by gamma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df = sym('df',[1 6]);
length_of_CVar = length(CVar1);
df_by_CVar = sym('df_by_CVar', [2, length_of_CVar]);
for j=1:length_of_CVar
    df_by_CVar(1,j) = complexdiff3(f, CVar1(j), 0);
    df_by_CVar(2,j) = complexdiff3(f, CVar1(j), 1);
end

tempRow = sym([0, 0, 0, 0, 0]);
for j=1:length_of_CVar
    if isreal(CVar1(j))==1
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*dCVar1(1,j);
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*dCVar1(2,j);
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*dCVar1(3,j);
        tempRow(4) = tempRow(4) + df_by_CVar(1,j)*dCVar1(4,j);
        tempRow(5) = tempRow(5) + df_by_CVar(1,j)*dCVar1(5,j);
    else
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*dCVar1(1,j)...
            + df_by_CVar(2,j)*conj(dCVar1(2,j));
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*dCVar1(2,j)...
            + df_by_CVar(2,j)*conj(dCVar1(1,j));
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*dCVar1(3,j)...
            + df_by_CVar(2,j)*conj(dCVar1(3,j));
        tempRow(4) = tempRow(4) + df_by_CVar(1,j)*dCVar1(4,j)...
            + df_by_CVar(2,j)*conj(dCVar1(5,j));
        tempRow(5) = tempRow(5) + df_by_CVar(1,j)*dCVar1(5,j)...
            + df_by_CVar(2,j)*conj(dCVar1(4,j));
    end
end

df(1) = tempRow(1);
df(2) = tempRow(2);
df(3) = tempRow(3);
df(4) = tempRow(4);
df(5) = tempRow(5);
df(6) = diff(f, gamma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%