syms p y
syms f df_p df_conjp df_y %df_conjp is real
syms d2f_pp d2f_pconjp d2f_conjpconjp d2f_py d2f_conjpy d2f_yy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
realVariable = [y, df_conjp, d2f_conjpy];
CVar = [p, y, f, df_p, df_conjp, df_y];
% derivativeDict
derivativeDict.p = [1; 0; 0];
derivativeDict.y = [0; 0; 1];
derivativeDict.f = [df_p; df_conjp; df_y];
derivativeDict.df_p = [d2f_pp; d2f_pconjp; d2f_py];
derivativeDict.df_conjp = [d2f_pconjp; d2f_conjpconjp; d2f_conjpy];
derivativeDict.df_y = [d2f_py; d2f_conjpy; d2f_yy];

MVar = [p, y, f, df_p, df_conjp, df_y, ...
    d2f_pp, d2f_pconjp, d2f_conjpconjp, d2f_py, d2f_conjpy, d2f_yy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% df = df_example_pconjpy(f,CVar,derivativeDict)
d2q_py = df_conjp*df_y - df_p*conj(df_y);
d2q_pconjp = 1/2*(df_conjp^2 -df_p*conj(df_p));
d2q_pp = -d2q_py^2/(4*d2q_pconjp);
d2q_yy = (d2q_py*conj(d2q_py)-4*d2q_pconjp^2)/(2*d2q_pconjp);
d2q_pp = complex_simple3(d2q_pp, MVar);
d2q_yy = complex_simple3(d2q_yy, MVar);

% find d3q_ppy and d3q_pyp
d3q_ppVec = df_example_pconjpy(d2q_pp, CVar, derivativeDict);
d3q_pyVec = df_example_pconjpy(d2q_py, CVar, derivativeDict);
d3q_ppy = d3q_ppVec(3);
d3q_ppy = complex_simple3(d3q_ppy, MVar);
d3q_pyp = d3q_pyVec(1);
d3q_pyp = complex_simple3(d3q_pyp, MVar);

part1 = (df_p*conj(df_y)-df_conjp*df_y)/(2*(df_p*conj(df_p)-df_conjp^2)^2);
term1 = (df_p*conj(df_y)-df_conjp*df_y)/(2*(df_p*conj(df_p)-df_conjp^2))...
    *(2*conj(d2f_yy)*df_p -2*d2f_yy*df_conjp +d2f_py*conj(df_y));
term2 = part1*(d2f_py*df_conjp*(conj(df_p)*df_y-df_conjp*conj(df_y))...
    +conj(d2f_py)*df_p*(df_y*df_conjp -conj(df_y)*df_p)...
    +2*d2f_conjpy*df_p*(df_conjp*conj(df_y)-conj(df_p)*df_y));

check_d3q_ppy = d3q_ppy - term1 - term2;
check_d3q_ppy = complex_simple3(check_d3q_ppy, MVar);


