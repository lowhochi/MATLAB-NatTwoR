syms u
mu1 =u^2-1; 
mu2= 2*u; 
mu3= i*(u^2+1);
v1 = i*(mu2*conj(mu3)-mu3*conj(mu2));
v2 = i*(mu3*conj(mu1)-mu1*conj(mu3));
v3 = i*(mu1*conj(mu2)-mu2*conj(mu1));
v1 = complex_simple3(v1,[u]);
v2 = complex_simple3(v2,[u]);
v3 = complex_simple3(v3,[u]);
norm_of_v = sqrt(v1*v1+v2*v2+v3*v3);
norm_of_v = complex_simple3(norm_of_v,[u]);
v1normv = v1/norm_of_v; %T1
v2normv = v2/norm_of_v; %T2
v3normv = v3/norm_of_v; %T3
v1normv = complex_simple3(v1normv, [u]);
v2normv = complex_simple3(v2normv, [u]);
v3normv = complex_simple3(v3normv, [u]);

Y= 1+u*conj(u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cf_g112_g221 = -1/Y^2*(conj(mu2)*v1normv+conj(mu1)*v2normv)*mu3*v3normv...
%     + 1/Y^4*mu1*mu2*conj(mu3)^2;
% cf_g112_g221  = expand(cf_g112_g221);
% cf_g112_g221 = complex_simple3(cf_g112_g221,[u]);
% 
% cf_g112_g312 = 1/Y^2*mu1*conj(mu3)...
%     + 1/Y^2*(conj(mu3)*v1normv+conj(mu1)*v3normv)*mu3*v3normv...
%     - 1/Y^4*mu1*mu3*conj(mu3)^2;
% cf_g112_g312  = expand(cf_g112_g312);
% cf_g112_g312 = complex_simple3(cf_g112_g312,[u]);
% 
% cf_g221_g312 = 1/Y^2*mu3*conj(mu2)...
%     - 1/Y^2*(conj(mu3)*v2normv+conj(mu2)*v3normv)*mu3*v3normv...
%     + 1/Y^4*mu2*mu3*conj(mu3)^2;
% cf_g221_g312  = expand(cf_g221_g312);
% cf_g221_g312 = complex_simple3(cf_g221_g312,[u]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cf_g112_g123 = 1/Y^2*mu1*conj(mu3)...
    +1/Y^2*(mu3*conj(mu1)*v1normv^2 +mu1*conj(mu1)*v3normv*v1normv)...
    -1/Y^4*mu1^2*conj(mu1)*conj(mu3);
cf_g112_g123  = expand(cf_g112_g123);
cf_g112_g123 = complex_simple3(cf_g112_g123,[u]);

cf_g112_g223 = -1/Y^2*conj(mu2)*mu3...
    + 1/Y^2*(mu3*conj(mu2)*v1normv^2 + mu1*conj(mu1)*v3normv*v2normv)...
    - 1/Y^4*mu2*mu1*conj(mu1)*conj(mu3);
cf_g112_g223 = expand(cf_g112_g223);
cf_g112_g223 = complex_simple3(cf_g112_g223,[u]);

cf_g112_g332 = -1/Y^2*(mu3*conj(mu3)*v1normv^2 + mu1*conj(mu1)*v3normv^2)...
    + 1/Y^4*mu1*mu3*conj(mu1)*conj(mu3);
cf_g112_g332 = expand(cf_g112_g332);
cf_g112_g332 = complex_simple3(cf_g112_g332,[u]);

cf_g221_g123 = -1/Y^2*(mu3*conj(mu1)*v1normv*v2normv+mu1*conj(mu2)*v1normv*v3normv)...
    +1/Y^4*mu1*mu2*conj(mu1)*conj(mu3);
cf_g221_g123 = expand(cf_g221_g123);
cf_g221_g123 = complex_simple3(cf_g221_g123,[u]);

cf_g221_g223 = -1/Y^2*(mu3*conj(mu2)*v1normv*v2normv+mu1*conj(mu2)*v3normv*v2normv)...
    + 1/Y^4*mu2^2*conj(mu1)*conj(mu3);
cf_g221_g223 = expand(cf_g221_g223);
cf_g221_g223 = complex_simple3(cf_g221_g223,[u]);

cf_g221_g332 = 1/Y^2*mu2*conj(mu1)...
    +1/Y^2*(mu3*conj(mu3)*v1normv*v2normv+mu1*conj(mu2)*v3normv*v3normv)...
    -1/Y^4*mu2*mu3*conj(mu1)*conj(mu3);
cf_g221_g332 = expand(cf_g221_g332);
cf_g221_g332 = complex_simple3(cf_g221_g332,[u]);

cf_g312_g123 = 1/Y^2*(mu3*conj(mu1)*v1normv*v3normv +mu1*conj(mu3)*v3normv*v1normv)...
    -1/Y^4*mu1*mu3*conj(mu1)*conj(mu3);
cf_g312_g123 = expand(cf_g312_g123);
cf_g312_g123 = complex_simple3(cf_g312_g123,[u]);

cf_g312_g223 = 1/Y^2*(mu3*conj(mu2)*v1normv*v3normv +mu1*conj(mu3)*v3normv*v2normv)...
    -1/Y^4*mu2*mu3*conj(mu1)*conj(mu3);
cf_g312_g223 = expand(cf_g312_g223);
cf_g312_g223 = complex_simple3(cf_g312_g223,[u]);

cf_g312_g332 = 1/Y^2*mu1*conj(mu3)...
    -1/Y^2*(mu3*conj(mu3)*v1normv*v3normv+ mu1*conj(mu3)*v3normv^2)...
    +1/Y^4*mu3^2*conj(mu1)*conj(mu3);
cf_g312_g332 = expand(cf_g312_g332);
cf_g312_g332 = complex_simple3(cf_g312_g332,[u]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%