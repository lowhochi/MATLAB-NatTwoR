function AInv = myInverse(A)
% A is an nxn square matrix
n = length(A);
M = sym('M',[n,n]); % M = Cofactor matrix of A
for j=1:n
    for k=1:n
        minorJK = find_minor(A,j,k);
        detOfMinorJK = mydet_fun(minorJK,5);
        temp = j+k;
        M(j,k) = ((-1)^temp)*detOfMinorJK; 
        clear minor_jk
        clear temp
        clear detOfMinorJK;
    end
end
detA = mydet_fun(A,n);
AInv = (1/detA)*transpose(M);