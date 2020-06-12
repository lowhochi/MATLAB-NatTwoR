function df = df_NatTwo_MuSet_CR_rmW(f, CVar, dCVar)
% Output: df = [df/dMu, df/dconjMu, df/dvnormv, df/du, df/dconju]
% CVar may contain both real and complex variables.
% dCVar(1,:) = differentiate CVar(:) by Mu;
% dCVar(2,:) = differentiate CVar(:) by conjMu;
% dCVar(3,:) = differentiate CVar(:) by vnormv;
% dCVar(4,:) = differentiate CVar(:) by u;
% dCVar(5,:) = differentiate CVar(:) by conju;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df = sym('df',[1 5]);
length_of_CVar = length(CVar);
df_by_CVar = sym('df_by_CVar', [2, length_of_CVar]);
for j=1:length_of_CVar
    df_by_CVar(1,j) = complexdiff3(f, CVar(j), 0);
    df_by_CVar(2,j) = complexdiff3(f, CVar(j), 1);
end

tempRow = sym([0, 0, 0, 0, 0]);
for j=1:length_of_CVar
    
    if isreal(CVar(j))==1
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*dCVar(1,j);
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*dCVar(2,j);
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*dCVar(3,j);
        tempRow(4) = tempRow(4) + df_by_CVar(1,j)*dCVar(4,j);
        tempRow(5) = tempRow(5) + df_by_CVar(1,j)*dCVar(5,j);
    else
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*dCVar(1,j)...
            + df_by_CVar(2,j)*conj(dCVar(2,j));
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*dCVar(2,j)...
            + df_by_CVar(2,j)*conj(dCVar(1,j));
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*dCVar(3,j)...
            + df_by_CVar(2,j)*conj(dCVar(3,j));
        tempRow(4) = tempRow(4) + df_by_CVar(1,j)*dCVar(4,j)...
            + df_by_CVar(2,j)*conj(dCVar(5,j));
        tempRow(5) = tempRow(5) + df_by_CVar(1,j)*dCVar(5,j)...
            + df_by_CVar(2,j)*conj(dCVar(4,j));
    end
end

df(1) = tempRow(1);
df(2) = tempRow(2);
df(3) = tempRow(3);
df(4) = tempRow(4);
df(5) = tempRow(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%