function df = df_Gijk_CR_rmW(f,GVar,GijkDict)
% Input: f = f(Gij_k, dGij_k_byl);
% Output: df = [df_by1, df_by2, df_by3];

% Don't use real variables as inpt!
df = sym('df',[1,3]);
length_of_GVar = length(GVar);
df_by1 = 0;
df_by2 = 0;
df_by3 = 0; 

df_by_GVar = sym('df',[1, length_of_GVar]);
for j=1:length_of_GVar
    df_by_GVar(j) = diff(f, GVar(j));
end

for j=1:length_of_GVar
    myChar = char(GVar(j));
    column = GijkDict.(myChar);
        df_by1 = df_by1 + df_by_GVar(j)*column(1);
        df_by2 = df_by2 + df_by_GVar(j)*column(2);
        df_by3 = df_by3 + df_by_GVar(j)*column(3);  
end

df(1) = df_by1;
df(2) = df_by2;
df(3) = df_by3;