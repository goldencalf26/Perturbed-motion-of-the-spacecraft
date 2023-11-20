clc; clear all; close all;
% Initial conditions
mu = 398600; % gravity coefficient
eps = 2.634e+10;% Earth compression constant
J = 1082e-6; % The second zonal harmonic
i_0 = deg2rad(75); % inclination
C_x = 2.2; % Coefficient of aerodynamic drag
unit = 3; % number of cubesat units
m = 3; % cubesat weight, kg
a = 0.1; % side length, m
S_mid = a^2;
sigma = sigma(C_x,S_mid,m); % ballistic coefficient
H_p0 = 390;
H_a0 = 410;
ro = density(H_p0,H_a0);
Omega_0 = deg2rad(10);
omega_0 = deg2rad(0);
thetta_0 = deg2rad(0);
R_e = 6371; % the radius of the Earth
Start_Conditions = StartCondition(H_p0, H_a0, thetta_0, R_e, omega_0, Omega_0, i_0,mu); % [x0,y0,z0,vx0,vy0,vz0]
t_interval = 1:1:60000;
% [t,X] = ode45(@equationofmotion,t_interval,Start_Conditions,odeset('RelTol',1e-8,'AbsTol',1e-10));
[t, X] = ode45(@(t, X) equationofmotion(t, X, sigma, ro, mu, J), t_interval, Start_Conditions, odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
r_mod = sqrt(X(:,1).^2 + X(:,2).^2 +  X(:,3).^2); % modulus of the radius of the vector
V_mod = sqrt(X(:,4).^2 + X(:,5).^2 +  X(:,6).^2); % modulus of the velocity

figure
plot(t_interval', r_mod,'red');
xlabel('t, s');
ylabel('r, km');
title('Modulus of the radius');
xlim([0 60000])
hold off
grid on

[Cxx, Cyy, Czz] = CoffC(X);

C = sqrt(Cxx.^2 + Cyy.^2 + Czz.^2);
Ci = Czz./C;

Om = 180 - mod(atan2d(-Cxx, Cyy), 360); % longitude of the ascending node
inc = mod(acosd(Czz./C), 180); % inclination
inc_w = acos(Czz./C);
u_w = asin(X(:,3)./(r_mod.*sin(inc_w)));

% projections of the disturbing acceleration on the radial, transversal and normal direction
[S_w, T_w, W_w] = STW(eps, r_mod, inc_w, u_w);

p_w = C.^2./mu;
nu_w = V_mod.*r_mod./mu;
tetta_wB = acos(sqrt(nu_w.*mu.*r_mod./C.^2));
e_w = sqrt(1 + nu_w.*(nu_w-2).*cos(tetta_wB).^2);
tetta_wm = atan(nu_w.*sin(tetta_wB).*cos(tetta_wB)./(nu_w.*cos(tetta_wB).^2 - 1));
w_w = sqrt(p_w./mu).*(-S_w.*cos(tetta_wm)./e_w + T_w.*(1 + r_mod./p_w).*sin(tetta_wm)./e_w - W_w.*r_mod./p_w.*cot(inc_w).*sin(u_w));...
    ...changing the argument of the pericenter at every second

% Plots
figure;
plot(t_interval, abs(Om),'red');
xlabel('Time, s');
ylabel('Longitude of the ascending node, deg');
title('Longitude of the ascending node');
ylim([0 180])
hold off
grid on

figure;
plot(t_interval, inc,'red');
xlabel('Time, s');
ylabel('Inclination, deg');
title('Inclination');
ylim([74 76])
hold off
grid on

w1 = omega_0*ones(size(t_interval));
w_list = w1' + w_w;
figure;
plot(t_interval, w_list,'red');
xlabel('Time, s');
ylabel('The pericenter argument, deg');
title('The pericenter argument');
hold off
grid on