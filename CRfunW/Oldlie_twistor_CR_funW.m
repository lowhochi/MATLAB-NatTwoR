function lie_funW = lie_twistor_CR_funW(T,u,w,dw,Mdw)
% M = [x, y, z, u]; 
% Mw= [u,w];
% dw = [dw_x, dw_y, dw_z, dw_u];
% M_dw = cat(2, Mw, dw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lie brackets without T. 
Dw_by_conj_mu=(conj(u)^2-1)*dw(1)+2*conj(u)*dw(2)-i*(conj(u)^2+1)*dw(3);
lie_X1_conj_X1 = [0; 0; 0; Dw_by_conj_mu; -conj(Dw_by_conj_mu)];
lie_X1_X2 = [0; 0; 0; 0; 0];
lie_X1_conj_X2 = [-2*conj(u); -2; 2*i*conj(u); 0; -conj(dw(4))];
lie_X2_conj_X1 = [2*u; 2; 2*i*u;dw(4);0]; 
lie_X2_conj_X2 = [0; 0; 0; 0; 0];
lie_conj_X1_conj_X2 = [0; 0; 0; 0; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DT_by_u and DT_by_conj_u.
% Lie brackets regarding T. 
DT_by_u1 = complexdiff3(T(1), u, 0);
DT_by_u2 = complexdiff3(T(2), u, 0);
DT_by_u3 = complexdiff3(T(3), u, 0);
DT_by_u1 = complex_simple3(DT_by_u1, Mdw);
DT_by_u2 = complex_simple3(DT_by_u2, Mdw);
DT_by_u3 = complex_simple3(DT_by_u3, Mdw);
%
DT_by_conj_u1 = complexdiff3(T(1), u, 1);
DT_by_conj_u2 = complexdiff3(T(2), u, 1);
DT_by_conj_u3 = complexdiff3(T(3), u, 1);
DT_by_conj_u1 = complex_simple3(DT_by_conj_u1, Mdw);
DT_by_conj_u2 = complex_simple3(DT_by_conj_u2, Mdw);
DT_by_conj_u3 = complex_simple3(DT_by_conj_u3, Mdw);
%
DT_by_u = [DT_by_u1; DT_by_u2; DT_by_u3; 0; 0];
DT_by_conj_u = [DT_by_conj_u1; DT_by_conj_u2; DT_by_conj_u3; 0; 0];
temp = T(1)*dw(1)+T(2)*dw(2)+T(3)*dw(3);
lie_X1_T = conj(w)*DT_by_conj_u - [0;0;0;0;conj(temp)];
lie_conj_X1_T = w*DT_by_u - [0;0;0;temp;0];
lie_X2_T = DT_by_u;
lie_conj_X2_T = DT_by_conj_u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lie_X1_X1= [0;0;0;0;0]; 
lie_conj_X1_conj_X1= [0;0;0;0;0];
lie_X2_X2= [0;0;0;0;0]; 
lie_conj_X2_conj_X2 = [0;0;0;0;0];
lie_T_T= [0;0;0;0;0];
% Output: lie_fun as a dictionary.
lie_funW = {lie_X1_X1, lie_X1_conj_X1, lie_X1_X2, lie_X1_conj_X2, lie_X1_T;
-lie_X1_conj_X1, lie_conj_X1_conj_X1, -lie_X2_conj_X1, lie_conj_X1_conj_X2, lie_conj_X1_T;
-lie_X1_X2, lie_X2_conj_X1, lie_X2_X2, lie_X2_conj_X2, lie_X2_T;
-lie_X1_conj_X2, -lie_conj_X1_conj_X2, -lie_X2_conj_X2, lie_conj_X2_conj_X2, lie_conj_X2_T;
-lie_X1_T, -lie_conj_X1_T, -lie_X2_T, -lie_conj_X2_T, lie_T_T};
