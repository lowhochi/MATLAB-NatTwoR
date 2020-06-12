function delta = checkArrayEqual(array1, array2)
n1 = length(array1);
n2 = length(array2);
count1 = 0;
count2 = 0;
if n1~=n2
    delta = 0;
else
    for j=1:n1
        temp = array1(j);
        if ismember(temp, array2)==1
            count1 = count1 + 1;
        end
    end
    
    for j=1:n2
        temp2 = array2(j);
        if ismember(temp2, array1)==1
            count2 = count2 + 1;
        end
    end
    
    if (n1==count1)&&(n2==count2)
        delta=1;
    else
        delta=0;
    end
end