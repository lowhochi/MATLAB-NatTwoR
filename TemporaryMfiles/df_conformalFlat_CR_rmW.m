function df = df_conformalFlat_CR_rmW(f, CVar, derivativeDict)
% output: df = [df/dx, df/dy, df/dz];
df = sym('df',[1,3]);
length_of_CVar = length(CVar);
df_by_CVar = sym('df_by_CVar', [2, length_of_CVar]);
for j=1:length_of_CVar
    df_by_CVar(1,j) = complexdiff3(f, CVar(j), 0);
    df_by_CVar(2,j) = complexdiff3(f, CVar(j), 1);
end
df(1) = 0;
df(2) = 0;
df(3) = 0;
for j=1:length_of_CVar
    myChar = char(CVar(j));
    columnTemp = derivativeDict.(myChar);  
    df(1) = df(1) +df_by_CVar(1,j)*columnTemp(1)...
        + df_by_CVar(2,j)*conj(columnTemp(1));
    df(2) = df(2) +df_by_CVar(1,j)*columnTemp(2)...
        + df_by_CVar(2,j)*conj(columnTemp(2));
    df(3) = df(3) +df_by_CVar(1,j)*columnTemp(3)...
        + df_by_CVar(2,j)*conj(columnTemp(3));
end