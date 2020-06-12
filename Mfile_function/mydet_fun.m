% Find the determinant of an nxn matrix A. 
function mydet = mydet_fun(A,n)

if n==1
    mydet = A(1);
    return;
elseif n==2
    mydet = A(1,1)*A(2,2)-A(1,2)*A(2,1);
    return;
end

temp_sum = 0;
for j=1:n
    Mj =find_minor(A,1,j);
    temp_sum = temp_sum + (-1)^(1+j)*A(1,j)*mydet_fun(Mj,n-1);
end
mydet = temp_sum;

% First row expansion of A.
% det(A)=a11*M11-a12*M12+...+(-1)^(n+1)*a1n*M1n




