function Mat2 = myRowOperation(Mat,rowChange,factorA,rowOp,factorB,MVarSet)
% In the matrix Mat, 
% Row operation on rowChange
% rowChange' = factorA*rowChange + factorB*rowOp
size_of_Mat = size(Mat);
N1 = size_of_Mat(1);
N2 = size_of_Mat(2);
Mat2 = sym('Mat2',[N1,N2]);
for j=1:N1
    for k=1:N2
        if (j==rowChange)
            Mat2(j,k) = factorA*Mat(j,k) +factorB*Mat(rowOp,k);
            Mat2(j,k) = complex_simple3(Mat2(j,k),MVarSet);
        else
            Mat2(j,k) = Mat(j,k);
        end
    end
end