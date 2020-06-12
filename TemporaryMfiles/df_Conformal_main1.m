function df = df_Conformal_main1(f,CVar,derivativeDict)
df = sym('df',[1 5]);
length_of_CVar = length(CVar);

df_by_CVar = sym('df_by_CVar', [2, length_of_CVar]);
for j=1:length_of_CVar
    df_by_CVar(1,j) = complexdiff3(f, CVar(j), 0);
    df_by_CVar(2,j) = complexdiff3(f, CVar(j), 1);
end

tempRow = sym([0, 0, 0, 0, 0]);
for j=1:length_of_CVar
    myChar = char(CVar(j));
    column = derivativeDict.(myChar); % Replace dCVar in future.
    if isreal(CVar(j))==1
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*column(1); %by Mu
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*column(2); %by conjMu
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*column(3); %by vnormv
        tempRow(4) = tempRow(4) + df_by_CVar(1,j)*column(4); %by u
        tempRow(5) = tempRow(5) + df_by_CVar(1,j)*column(5); %by conj(u)
    else
        tempRow(1) = tempRow(1) + df_by_CVar(1,j)*column(1)...
            + df_by_CVar(2,j)*conj(column(2)); %by Mu
        tempRow(2) = tempRow(2) + df_by_CVar(1,j)*column(2)...
            + df_by_CVar(2,j)*conj(column(1)); %by conjMu
        tempRow(3) = tempRow(3) + df_by_CVar(1,j)*column(3)...
            + df_by_CVar(2,j)*conj(column(3)); %by vnormv 
        tempRow(4) = tempRow(4) + df_by_CVar(1,j)*column(4)...
            + df_by_CVar(2,j)*conj(column(5)); %by u 
        tempRow(5) = tempRow(5) + df_by_CVar(1,j)*column(5)...
            + df_by_CVar(2,j)*conj(column(4)); %by conj(u)
    end
end

df(1) = tempRow(1);
df(2) = tempRow(2);
df(3) = tempRow(3);
df(4) = tempRow(4);
df(5) = tempRow(5);