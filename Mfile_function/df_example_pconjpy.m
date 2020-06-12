function df = df_example_pconjpy(f,CVar,derivativeDict)
% input f = f(p,conj(p),y)
% output [df/dp, df/dconjp, df/dy];
df = sym('df',[1 3]);
length_of_CVar = length(CVar);

df_by_CVar = sym('df_by_CVar', [2, length_of_CVar]);
for j=1:length_of_CVar
    df_by_CVar(1,j) = complexdiff3(f, CVar(j), 0);
    df_by_CVar(2,j) = complexdiff3(f, CVar(j), 1);
end

tempRow = sym([0, 0, 0]);
for j=1:length_of_CVar
    myChar = char(CVar(j));
    column = derivativeDict.(myChar); % Replace dCVar in future.
    if isreal(CVar(j))==1
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*column(1); %by p
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*column(2); %by conjp
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*column(3); %by y
    else
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*column(1)...
            + df_by_CVar(2,j)*conj(column(2)); %by p
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*column(2)...
            + df_by_CVar(2,j)*conj(column(1)); %by conjp
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*column(3)...
            + df_by_CVar(2,j)*conj(column(3)); %by y 
    end
end
df(1) = tempRow(1);
df(2) = tempRow(2);
df(3) = tempRow(3);