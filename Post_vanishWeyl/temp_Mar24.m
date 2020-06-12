% load('DataWeyl7_Part25_NatTwoR_CR_rmW.mat');
% variable_NatTwoR_W725_CR_rmW
% clearvars myChar myCharR subChar wChar
% % % myWorkSpace = who; %[1100,1]
% % % d2GList = [];
% % % for j=1:length(myWorkSpace)
% % %     tempCell = myWorkSpace(j);
% % %     tempChar = tempCell{1};
% % %     tempString = convertCharsToStrings(tempChar);
% % %     if strfind(tempString, 'd2G')==1
% % %         d2GList = [d2GList; tempString];
% % %     end
% % % end
% % % % d2GList: string array, 292x1;
% % % d2Gset = [];
% % % for j=1:length(d2GList)
% % %     myString = d2GList(j);
% % %     if (strfind(myString,'d2G')==1)
% % %         if (strfind(myString,'by')==9)
% % %             if (strfind(myString,'Sub')==13)
% % %                 continue
% % %             else
% % %                 disp(myString);
% % %                 eval(['d2Gset=[d2Gset,',char(myString),'];']);
% % %             end
% % %         end
% % %     end
% % % end
% % index of covPSety
% indexArray = ["112", "113", "122", "123", "133", ...
%     "211", "212", "213", "223", "233", "311", "312", ...
%     "313", "322", "323"];
% % covPSet = [covP112, covP113, covP122, covP123, covP133,...
% %     covP211, covP212, covP213, covP223, covP233,...
% %     covP311, covP312, covP313, covP322, covP323];
% for j=1:15
%     char01 = char(indexArray(j));
%     m = str2num(char01(1));
%     n = str2num(char01(2));
%     k = str2num(char01(3));
%     [termVec, gVec] = coeffs(covP(m,n,k),d2Gset);
%     disp(gVec(end));
%     eval(['Aterm',char01,'=termVec(end);'])
% end

% %% 
syms u y C112 C113 C221 C321 C123
w1212 = 4*i*y^2*(conj(u)^2-u^2)*C112 ...
    + 4*y^2*(u+conj(u))*(1-u*conj(u))*C113 ...
    + 4*i*y^2*(u-conj(u))*(1-u*conj(u))*C221 ...
    - 2*y^2*((u-conj(u))^2+(1-u*conj(u))^2)*C321 ...
    - 2*y^2*((u+conj(u))^2-(1-u*conj(u))^2)*C123;
    
[cf1212, cVec] = coeffs(w1212, [C112, C113, C221, C321, C123]);
for j=1:length(cf1212)
    temp = expand(cf1212(j)/(4*y^2));
    disp(temp);
end








