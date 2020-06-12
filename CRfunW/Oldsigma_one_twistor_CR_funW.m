function sigma0 = sigma_one_twistor_CR_funW(Theta, ricci_curv, h, omega, ...
    conj_omega,Md2w)
% (1) Part 1: \sum_{n=1}^2 omega_n^n.
% (2) Part 2: -(i/2)*(h^{\alpha \beta}*dh_{\alpha \beta}).
% (3) Part 3: -(1/12)*rho*alpha0.
theta = Theta(5,:);
rho = ricci_curv{2};
sigma0_part1 = i*omega{1,1}+i*omega{2,2};
sigma0_part3 = (-1/12)*rho*theta;
% The differential of h.
% dh_{m \on}= \sum_{k=1}^2 h_{m \ok}*conj(\omega_n^k)+ h_{k\on}*\omega_m^k. 
h_inv = inv(h);
dh11 = zeros(1,5);
dh12 = zeros(1,5);
dh21 = zeros(1,5);
dh22 = zeros(1,5);
for k=1:2
    dh11 = dh11 + h(1,k)*conj_omega{1,k} + h(k,1)*omega{1,k};
    dh12 = dh12 + h(1,k)*conj_omega{2,k} + h(k,2)*omega{1,k};
    dh21 = dh21 + h(2,k)*conj_omega{1,k} + h(k,1)*omega{2,k};
    dh22 = dh22 + h(2,k)*conj_omega{2,k} + h(k,2)*omega{2,k};
end
Dh = {dh11, dh12; dh21, dh22};
sigma0_part2 = (-i/2)*(h_inv(1,1)*Dh{1,1} + h_inv(2,1)*Dh{1,2} ...
    + h_inv(1,2)*Dh{2,1} + h_inv(2,2)*Dh{2,2}); % Change in h here.
%
sigma0_part2_1 = complex_simple3(sigma0_part2(1), Md2w);
sigma0_part2_2 = complex_simple3(sigma0_part2(2), Md2w);
sigma0_part2_3 = complex_simple3(sigma0_part2(3), Md2w);
sigma0_part2_4 = complex_simple3(sigma0_part2(4), Md2w);
sigma0_part2_5 = complex_simple3(sigma0_part2(5), Md2w);
sigma0_part2_check = [sigma0_part2_1, sigma0_part2_2, sigma0_part2_3, ...
    sigma0_part2_4, sigma0_part2_5];
% disp(sigma0_part2_check);
%
sigma0 = sigma0_part1 + sigma0_part2 + sigma0_part3; 




