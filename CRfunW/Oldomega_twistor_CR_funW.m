function omega = omega_twistor_CR_funW(Theta, Gamma)
omega0 = sym('omega',[2 2 5]);
% omega(m,n,:} = omega_m^n as a 1-form on (N,CR_W) 
% omega_m^n = Gamma_{km}^n*theta^k + Gamma_(\ok m}^n*theta^{\ok}
%               + Gamma_{0m}^n*alpha0.
Theta_holo = [Theta(1,:); Theta(3,:)];
Theta_antiholo = [Theta(2,:); Theta(4,:)];
for m=1:2
    for n=1:2
        part1 = zeros([1,5]);
        part2 = zeros([1,5]);
        for k=1:2
            part1 = part1 + Gamma.holo(k,m,n)*Theta_holo(k,:);
            part2 = part2 + Gamma.antiholo(k,m,n)*Theta_antiholo(k,:);
        end
        temp_1_form = part1 + part2 + Gamma.T(m,n)*Theta(5,:);
        omega0(m,n,:) = temp_1_form;
        clearvars part1 part2
        clear temp_1_form
    end
end

omega11 = sym('omega11',[1,5]);
omega12 = sym('omega12',[1,5]);
omega21 = sym('omega21',[1,5]);
omega22 = sym('omega22',[1,5]);
for k=1:5
    omega11(k) = omega0(1,1,k);
    omega12(k) = omega0(1,2,k);
    omega21(k) = omega0(2,1,k);
    omega22(k) = omega0(2,2,k);
end
% Every omega11 & etc is a 1x5 row vector. 
omega = {omega11, omega12; omega21, omega22};
                
            
            