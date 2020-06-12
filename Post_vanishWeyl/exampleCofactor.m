function [d2Qfun, Cmat]= exampleCofactor(Qfun, varSet)
% Qfun depends on p, conjp, y, or other
% varSet = [p, y];
dQ_p = complexdiff3(Qfun, varSet(1), 0);
dQ_y = diff(Qfun, varSet(2));
d2Q_pp = complexdiff3(dQ_p, varSet(1), 0);
d2Q_pconjp = complexdiff3(dQ_p, varSet(1), 1);
d2Q_py = diff(dQ_p, varSet(2));
d2Q_yy = diff(dQ_y, varSet(2));

d2Qfun = [d2Q_pp, d2Q_pconjp, d2Q_py;
    d2Q_pconjp, conj(d2Q_pp), conj(d2Q_py);
    d2Q_py, conj(d2Q_py), d2Q_yy];

Cmat = [2*conj(d2Q_pp), -2*d2Q_pconjp-d2Q_yy, conj(d2Q_py);
    -2*d2Q_pconjp-d2Q_yy, 2*d2Q_pp, d2Q_py;
    conj(d2Q_py), d2Q_py, -2*d2Q_pconjp];