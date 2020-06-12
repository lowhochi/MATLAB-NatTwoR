function Cmn = myCofactor(A, ii, jj)

N = length(A);
A0 = sym('A0',[N-1,N-1]);

rowCount = 1;
columnCount =1;

for m = 1:N
    if m==ii
        continue
    end
    for n=1:N
        if (n==jj)
            continue
        end
        A0(rowCount,columnCount) = A(m,n);
        columnCount = columnCount+1;
    end   
    rowCount = rowCount + 1;
    columnCount = 1;
end

Cmn = (-1)^(ii+jj)*det(A0);