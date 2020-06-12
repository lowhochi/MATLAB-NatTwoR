load('Data_solveD2q_detPhiZero.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[d2Qfun, CoQ]=exampleCofactor(Qfun, varSet, MVar)
syms a0 real
syms N integer
Qfun03 = -p^2-conj(p)^2-2*y^2;
[d2Q03, Cmat03] = exampleCofactor(Qfun03, [p,y]);
totalQ03 = sym('totalQ03',[9,6]);
% d2qVec =[d2q_pp, d2q_py, d2q_pconjp, d2q_yy];
d2QVec03 = [d2Q03(1,1), d2Q03(1,3), d2Q03(1,2), d2Q03(3,3)];
det_of_Cmat03 = det(Cmat03);
det_of_Cmat03 = complex_simple3(det_of_Cmat03, MVar);
for j=1:9
   for k=1:6
       temp = totalMat(j,k);
       temp = subs(temp, d2qVec, d2QVec03);
       totalQ03(j,k) = complex_simple3(temp, MVar);
   end
end
clearvars j k temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qfun01 = a0*p^N + a0*conj(p)^N;
% [d2Q01, CoQ01] = exampleCofactor(Qfun01, [p,y], MVar);
% totalQ01 = sym('totalQ01',[9,6]);
% % d2qVec =[d2q_pp, d2q_py, d2q_pconjp, d2q_yy];
% d2QVec01 = [d2Q01(1,1), d2Q01(1,3), d2Q01(1,2), d2Q01(3,3)];
% for j=1:9
%    for k=1:6
%        temp = totalMat(j,k);
%        temp = subs(temp, d2qVec, d2QVec01);
%        totalQ01(j,k) = complex_simple3(temp, MVar);
%    end
% end
% clearvars j k temp
% subTotalQ01 = totalQ01(:,1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qfun02 = a0*p*y +a0*conj(p)*y;
% [d2Q02, CoQ02] = exampleCofactor(Qfun02, [p,y], MVar);
% totalQ02 = sym('totalQ02',[9,6]);
% det_of_d2Q02 = det(d2Q02);
% det_of_d2Q02 = complex_simple3(det_of_d2Q02, MVar);
% % d2qVec =[d2q_pp, d2q_py, d2q_pconjp, d2q_yy];
% d2QVec02 = [d2Q02(1,1), d2Q02(1,3), d2Q02(1,2), d2Q02(3,3)];
% for j=1:9
%    for k=1:6
%        temp = totalMat(j,k);
%        temp = subs(temp, d2qVec, d2QVec02);
%        totalQ02(j,k) = complex_simple3(temp, MVar);
%    end
% end
% clearvars j k temp
% subTotalQ02 = [totalQ02(:,1),totalQ02(:,2),totalQ02(:,5)];
% Cmat02 = sym('Cmat02',[3,3]);
% for j=1:3
%     for k=1:3
%         Cmat02(j,k) = subs(Cmat(j,k), d2qVec, d2QVec02);
%         Cmat02(j,k) = complex_simple3(Cmat02(j,k), MVar);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
