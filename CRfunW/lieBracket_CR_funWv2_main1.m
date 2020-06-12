function Z = lieBracket_CR_funWv2_main1(X,Y,CVar,derivativeDict)
% Given X = a1*d/dx + a2*d/dy + a3*d/dz + a4*d/du +a5*d/dconj(u),
% Y = b1*d/dx + b2*d/dy + b3*d/dz + b4*d/du +b5*d/dconj(u),
% find [X,Y] in [d/dx; d/dy; d/dz; d/du; d/dconj(u)].
dX = sym('dY',[5 5]);
dY = sym('dY',[5 5]);
% dX(m,n) = dX(m) by x_n (or u, conj(u)).
for m=1:5
    dXm = df_main1_CR_funWv2(X(m),CVar,derivativeDict);
    dYm = df_main1_CR_funWv2(Y(m),CVar,derivativeDict);
    for n=1:5
        dX(m,n)=dXm(n);
        dY(m,n)=dYm(n);
    end
end

Z = sym('Z',[5,1]); %Z=[X,Y]
for jj=1:5
    temp1 = 0;
    for ii=1:5
        temp1 = temp1 + X(ii)*dY(jj,ii) - Y(ii)*dX(jj,ii);
    end
    Z(jj) = temp1;
end