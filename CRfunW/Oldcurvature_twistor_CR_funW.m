function R = curvature_twistor_CR_funW(g1,h,Gamma,holo_vector,antiholo_vector, ...
    lie_funW,T,RVar2,CVar2,dRVar2,dCVar2,MVar2)

% The curvature tensor R is in 2 x 2 x 2 x 2.
% R(m,n,k,l) = R_m^n_k_ol where R(X_k, conj_X_l)X_m = R_m^n_k_ol X_n. 
variable = {'holo_x','holo_y','holo_z','holo_u','holo_conj_u';
    'antiholo_x','antiholo_y','antiholo_z','antiholo_u','antiholo_conj_u';
    'T_x','T_y','T_z','T_u','T_conj_u'};
%
dGamma.holo_x = sym('dGamma1_%d_%d_%d', [2,2,2]);
dGamma.holo_y = sym('dGamma2_%d_%d_%d', [2,2,2]);
dGamma.holo_z = sym('dGamma3_%d_%d_%d', [2,2,2]);
dGamma.holo_u = sym('dGamma4_%d_%d_%d', [2,2,2]);
dGamma.holo_conj_u = sym('dGamma5_%d_%d_%d', [2,2,2]);
%
dGamma.antiholo_x = sym('dGamma7_%d_%d_%d', [2,2,2]);
dGamma.antiholo_y = sym('dGamma8_%d_%d_%d', [2,2,2]);
dGamma.antiholo_z = sym('dGamma9_%d_%d_%d', [2,2,2]);
dGamma.antiholo_u= sym('dGamma10_%d_%d_%d', [2,2,2]);
dGamma.antiholo_conj_u = sym('dGamma11_%d_%d_%d', [2,2,2]);
%
dGamma.T_x = sym('dGamma13_%d_%d', [2,2]);
dGamma.T_y = sym('dGamma14_%d_%d', [2,2]);
dGamma.T_z = sym('dGamma15_%d_%d', [2,2]);
dGamma.T_u = sym('dGamma16_%d_%d', [2,2]);
dGamma.T_conj_u = sym('dGamma17_%d_%d', [2,2]);
%
for m=1:2
    for n=1:2
        for k=1:2
            temp1_fun = Gamma.holo(m,n,k);
            temp2_fun = Gamma.antiholo(m,n,k);
            temp1_dGamma= df_twistor_CR_funW(temp1_fun,RVar2,CVar2,dRVar2,dCVar2,MVar2);
            temp2_dGamma= df_twistor_CR_funW(temp2_fun,RVar2,CVar2,dRVar2,dCVar2,MVar2);
            dGamma.holo_x(m,n,k) = temp1_dGamma(1);
            dGamma.holo_y(m,n,k) = temp1_dGamma(2);
            dGamma.holo_z(m,n,k) = temp1_dGamma(3);
            dGamma.holo_u(m,n,k)= temp1_dGamma(4);
            dGamma.holo_conj_u(m,n,k)= temp1_dGamma(5);
            %
            dGamma.antiholo_x(m,n,k) = temp2_dGamma(1);
            dGamma.antiholo_y(m,n,k) = temp2_dGamma(2);
            dGamma.antiholo_z(m,n,k) = temp2_dGamma(3);
            dGamma.antiholo_u(m,n,k)= temp2_dGamma(4);
            dGamma.antiholo_conj_u(m,n,k)= temp2_dGamma(5);
            %
            clear temp1_fun
            clear temp2_fun
            clear temp1_dGamma
            clear temp2_dGamma
        end
    end
end

for n=1:2
    for k=1:2
        temp3_fun = Gamma.T(n,k);
        temp3_dGamma= df_twistor_CR_funW(temp3_fun,...
            RVar2,CVar2,dRVar2,dCVar2,MVar2);
        dGamma.T_x(n,k)= temp3_dGamma(1);
        dGamma.T_y(n,k)= temp3_dGamma(2);
        dGamma.T_z(n,k)= temp3_dGamma(3);
        dGamma.T_u(n,k)= temp3_dGamma(4);
        dGamma.T_conj_u(n,k)= temp3_dGamma(5);
        clear temp3_fun
        clear temp3_dGamma
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(2) combined with Part(3):
R = sym('R0_%d_%d_%d_%d', [2,2,2,2]);
h_inv = inv(h);
for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                dGamma_by_X_k = 0;
                dGamma_by_X_ol = 0;
                Gamma_sum = 0;
                term_a = 0; %p3
                term_b = 0; %p3
                V_k = holo_vector(:,k); %p2, p3
                V_ll = antiholo_vector(:,ll); %p2, p3
                % Replace lie_bracket.
                % lie_vector = lie_bracket_twistor_CR_model(V_k,V_ll,M); %p3
                lie_vector = lie_funW{2*k-1, 2*ll};
                c = transpose(T)*g1*lie_vector; %p3
                term_c = Gamma.T(m,n)*c; %p3
                % Part(2) core:
                dG_ol_m_n=[dGamma.(variable{2,1})(ll,m,n), dGamma.(variable{2,2})(ll,m,n), ...
                    dGamma.(variable{2,3})(ll,m,n), dGamma.(variable{2,4})(ll,m,n), ...
                    dGamma.(variable{2,5})(ll,m,n)];
                dG_k_m_n=[dGamma.(variable{1,1})(k,m,n), dGamma.(variable{1,2})(k,m,n), ...
                    dGamma.(variable{1,3})(k,m,n), dGamma.(variable{1,4})(k,m,n), ...
                    dGamma.(variable{1,5})(k,m,n)];
                dGamma_by_X_k = dG_ol_m_n*V_k;
                dGamma_by_X_ol = dG_k_m_n*V_ll;
                for q=1:2
                    Gamma_sum = Gamma_sum + Gamma.antiholo(ll,m,q)*Gamma.holo(k,q,n) ...
                        -Gamma.holo(k,m,q)*Gamma.antiholo(ll,q,n);
                end
                % Part(2) core ends.
                % Part(3) core:
                for pp=1:2
                    a_p = 0;
                    b_p = 0;
                    for qq=1:2
                        V_q = holo_vector(:,qq);
                        W_q = antiholo_vector(:,qq);
                        a_p_0 = transpose(W_q)*g1*lie_vector;
                        b_p_0 = transpose(V_q)*g1*lie_vector;
                        a_p = a_p + h_inv(qq,pp)*a_p_0;
                        b_p = b_p + h_inv(pp,qq)*b_p_0;
                        clearvars V_q W_q a_p_0 b_p_0
                    end
                    term_a = term_a + a_p*Gamma.holo(pp,m,n);
                    term_b = term_b + b_p*Gamma.antiholo(pp,m,n);
                    clearvars a_p b_p
                end
                % Part(3) core ends. 
                % Define R(m,n,k,ll).
                R(m,n,k,ll) = dGamma_by_X_k - dGamma_by_X_ol + Gamma_sum ...
                    -term_a-term_b-term_c;
                % R(m,n,k,ll) = simplify(R(m,n,k,ll));
                % Clear variables. 
                clearvars Gamma_sum dGamm_by_X_ol dGamma_by_X_k
                clearvars dG_k_m_n dG_ol_m_n V_ll V_k
                clearvars term_a term_b c term_c lie_vector
            end
        end
    end
end