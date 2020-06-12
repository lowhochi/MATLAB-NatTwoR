% sigma is a 1-form in [\theta^1, \theta^\bar1, \theta^2, \theta^\bar2, \alpha, d\gamma];
% dSigma = exterior derivatives of sigma.
% When sigma = a_i*u_i*,
% dSigma = da_i\wedge u_i^* + a_i*du_i*;
function dSigma = df_flatM3_CR_rmW(sigma, dCoframe, CVar, derivativeDict, gamma)
w1 = sym('w1');
for kk=1:length(CVar)
    myChar = char(CVar(kk));
    if myChar=='w'
        w1=CVar(kk);
        break
    end
end

part1 = zeros(6,6);
% find da_i
for j=1:6
    aI = sigma(j);
    daI = df_main4_CR_funWv2(aI,CVar,derivativeDict,gamma);
    % in terms of [\theta^1, \theta^\bar1, \theta^2, \theta^\bar2, \alpha, d\gamma];
    daITwo = [daI(2)+conj(w1)*daI(5), daI(1)+w1*daI(4), daI(3), daI(4), daI(5), daI(6)];
% construct da_i\wedge u_i^*
    daI_Times_uI = sym('part1',[6 6]);
    uI_Times_daI = sym('part2',[6 6]);
    for m=1:6
        for n=1:6
            if m==j
                uI_Times_daI(m,n) = daITwo(n);
            else
                uI_Times_daI(m,n) = 0;
            end

            if n==j
                daI_Times_uI(m,n) = daITwo(m);
            else
                daI_Times_uI(m,n) = 0;
            end
        end
    end
    part1 = part1 + 1/2*(daI_Times_uI - uI_Times_daI);
end

part2 = zeros(6,6);
for j=1:6
    part2 = part2 + sigma(j)*dCoframe{1,j};
end
dSigma = part1 + part2;