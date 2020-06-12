% B is a matrix of size n by n.(n>1)
% 1<=r<=n and 1<=s<=n.
% Return the minor matrix of B at position (r,s).
function minor = find_minor(B,r,s)
size_of_B = size(B);
n = size_of_B(1);
minor = [];
% k is the column index of matrix B.
for k=1:n
    temp_column2 = [];
    if k~=s
        % j is the row index of matrix B.
        for j=1:n
            if j~=r
                temp_column2 = [temp_column2; B(j,k)];
            end
        end
        minor = [minor, temp_column2];
    end
end
