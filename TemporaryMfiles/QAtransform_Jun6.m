syms lambda phi s t 
syms u u2 % u2 is an alternate to uPrime
A = [lambda, -conj(phi);
    phi, conj(lambda)];
stVec = [s; t];
PrimeVec = A*stVec;
sPrime = PrimeVec(1); tPrime = PrimeVec(2);
sPrime = subs(sPrime,[s t],[u 1]);
tPrime = subs(tPrime,[s t],[u 1]);
uPrime = sPrime/tPrime;
uPrime = subs(uPrime,[s t],[u 1]);
%
qA = sym('qA',[3 3]);
qA(1,1) = 1/2*(lambda^2-phi^2+conj(lambda)^2-conj(phi)^2);
qA(1,2) = -lambda*conj(phi)-conj(lambda)*phi;
qA(1,3) = -i/2*(lambda^2-phi^2-conj(lambda)^2+conj(phi)^2);
qA(2,1) = lambda*phi + conj(lambda)*conj(phi);
qA(2,2) = lambda*conj(lambda) - phi*conj(phi);
qA(2,3) = -i*(lambda*phi - conj(lambda)*conj(phi));
qA(3,1) = i/2*(lambda^2+phi^2-conj(lambda)^2-conj(phi)^2);
qA(3,2) = i*(-lambda*conj(phi)+conj(lambda)*phi);
qA(3,3) = 1/2*(lambda^2+phi^2+conj(lambda)^2+conj(phi)^2);
%
muVec = [s^2-t^2; 2*s*t; i*(s^2+t^2)];
muPrimeVec = [sPrime^2- tPrime^2;
    2*sPrime*tPrime;
    i*(sPrime^2 + tPrime^2)]; % muPrimeVec = qA*muVec;
uPrimeVec = 1/tPrime^2*muPrimeVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test One
syms c0 c1 z0
syms k real
w = c0 + c1*u + k*u^2 - conj(c1)*u^3 + conj(c0)*u^4;
wPrime = conj(c0)*(1+u2*conj(z0))^4 ...
    - conj(c1)*(1+u2*conj(z0))^3*(z0-u2) ...
    + k*(1+u2*conj(z0))^2*(z0-u2)^2 ...
    + c1*(1+u2*conj(z0))*(z0-u2)^3 ...
    + c0*(z0-u2)^4;

[wPcoeff, wPterm] = coeffs(wPrime, [u2]);

