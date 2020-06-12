% twistor_Conformal_CR_main3.m
load('Data_Conformal_CR_Main2.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChernC = sym('ChernC',[2,2,2,2]);
RcThree = sym('RcThree',[2 2 2 2]);
for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                 RcThree(m,n,k,ll)=Rc(m,1,k,ll)*hC(1,n)+Rc(m,2,k,ll)*hC(2,n);
                 %
                 tempPart1 = RcThree(m,n,k,ll);
                 tempPart2 = (-1/4)*(hC(k,ll)*ricC(m,n) + hC(m,ll)*ricC(k,n)...
                    + ricC(k,ll)*hC(m,n) + ricC(m,ll)*hC(k,n));
                 tempPart3 = (rhoC/12)*(hC(k,ll)*hC(m,n) + hC(m,ll)*hC(k,n));
                 tempChern = tempPart1 + tempPart2 + tempPart3;
                 tempChern = subs(tempChern, conjRealVariable, realVariable);
                 ChernC(m,n,k,ll) = complex_simple3(tempChern,MVarRc);
            end
        end
    end
end
clearvars m n k ll tempPart1 tempPart2 tempPart3 tempChern
save('Data_Conformal_CR_Main3.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%