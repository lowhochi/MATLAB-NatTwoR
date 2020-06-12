function Gamma=christoffel_twistor_CR_funW(g1,h,holo_vector,antiholo_vector,T,lie_funW,M)
% \nabla_{X_m} X_n = \Gamma_{mn}^k X_k
% \nabla_{X_{\overline m}} X_n = Gamma_{{\overline m} n}^k X_k
% \nabla_T variable = {'holo', 'antiholo', 'T'};
% Gamma.holo(m,n,k) = \Gamma_{mn}^k
% Gamma.antiholo(m,n,k) = Gamma_{{\overline m} n}^k 
% Gamma.T(n,k) = Gamma_{0n}^k
variable = {'holo', 'antiholo', 'T'};
Gamma.holo = sym('Gamma_holo_%d_%d_%d',[2 2 2]);
Gamma.antiholo = sym('Gamma_antiholo_%d_%d_%d',[2 2 2]);
Gamma.T = sym('Gamma_T_%d_%d', [2 2]);
%
% The differential of h.
h_inv = inv(h);
dh.x = [complexdiff3(h(1,1),M(1),0), complexdiff3(h(1,2),M(1),0);
    complexdiff3(h(2,1),M(1),0), complexdiff3(h(2,2),M(1),0)];
dh.y = [complexdiff3(h(1,1),M(2),0), complexdiff3(h(1,2),M(2),0);
    complexdiff3(h(2,1),M(2),0), complexdiff3(h(2,2),M(2),0)];
dh.z = [complexdiff3(h(1,1),M(3),0), complexdiff3(h(1,2),M(3),0);
    complexdiff3(h(2,1),M(3),0), complexdiff3(h(2,2),M(3),0)];
dh.u = [complexdiff3(h(1,1),M(4),0), complexdiff3(h(1,2),M(4),0);
    complexdiff3(h(2,1),M(4),0), complexdiff3(h(2,2),M(4),0)];
dh.conj_u = [complexdiff3(h(1,1),M(4),1), complexdiff3(h(1,2),M(4),1);
    complexdiff3(h(2,1),M(4),1), complexdiff3(h(2,2),M(4),1)];
%
% dh_holo(m,n,l) = dh_{nl}(X_m).
dh_holo = sym('dh_%d_%d_%d', [2 2 2]);
for m=1:2
    for n=1:2
        for ll=1:2
            temp_vector = holo_vector(:,m);
            temp_dh = [dh.x(n,ll), dh.y(n,ll), dh.z(n,ll), ...
                dh.u(n,ll), dh.conj_u(n,ll)];
            dh_holo(m,n,ll) = temp_dh*temp_vector;
            clear temp_vector
            clear temp_dh
        end
    end
end
% 
% Define Gamma.holo(m,n,k). 
for m=1:2
    for n=1:2
        V_m = holo_vector(:,m);
        V_n = holo_vector(:,n);
        for k=1:2
            temp_part1 = 0;
            temp_part2 = 0;
            for ll=1:2
                % lie_vector = lie_bracket_twistor_CR_model(V_m, V_ll, M);
                lie_vector = lie_funW{2*m-1,2*ll};
                temp_part2_0 = transpose(V_n)*g1*lie_vector;
                temp_part2 = temp_part2 - h_inv(ll,k)*temp_part2_0;
                temp_part1 = temp_part1 + h_inv(ll,k)*dh_holo(m,n,ll);
                clear temp_part2_0
                clear lie_vector
            end
            Gamma.holo(m,n,k)= temp_part1 + temp_part2;
            % Gamma.holo(m,n,k)= simplify(Gamma.holo(m,n,k));
            clear temp_part1
            clear temp_part2
        end
        clear V_n
        clear V_m
    end
end

% Define Gamma.antiholo(m,n,k).
for m=1:2
    for n=1:2
        V_m = antiholo_vector(:,m);
        V_n = holo_vector(:,n);
        % lie_vector = lie_bracket_twistor_CR_model(V_m, V_n, M);
        lie_vector = lie_funW{2*m,2*n-1};
        for k=1:2
            temp_antiholo = 0;
            for ll=1:2
                V_ll = antiholo_vector(:,ll);
                temp_antiholo_0 = transpose(V_ll)*g1*lie_vector;
                temp_antiholo = temp_antiholo + h_inv(ll,k)*temp_antiholo_0;
                clear V_11
                clear temp_antiholo_0
            end
            Gamma.antiholo(m,n,k) = temp_antiholo;
            % Gamma.antiholo(m,n,k)= simplify(Gamma.antiholo(m,n,k));
            clear temp_antiholo
        end
        clear lie_vector
        clear V_n
        clear V_m
    end
end

% Define Gamma.T(n,k)
for n=1:2
    V_n = holo_vector(:,n);
    % lie_vector = lie_bracket_twistor_CR_model(T, V_n, M);
    lie_vector = lie_funW{5,2*n-1};
    for k=1:2    
        temp_T = 0;
        for ll=1:2
            V_ll = antiholo_vector(:,ll);
            temp_T_0 = transpose(V_ll)*g1*lie_vector;
            temp_T = temp_T + h_inv(ll,k)*temp_T_0;
            clear temp_T_0
            clear V_ll
        end
        Gamma.T(n,k) = temp_T;
        % Gamma.T(n,k) = simplify(Gamma.T(n,k));
        clear temp_T
    end
    clear V_n
    clear lie_vector
end