function SC = StartCondition(H_p0, H_a0, thetta_0, R_e, omega_0, Omega_0, i_0,mu)
r_p = R_e + H_p0;
r_a = R_e + H_a0;
a = (r_a + r_p)/2;
e = (r_a - r_p)/(r_a + r_p);
p = r_p*(1+e);
u = omega_0 + thetta_0;
r0 = p/(1 + e*cos(thetta_0));
x0 = r0*(cos(u)*cos(Omega_0) - sin(u)*cos(i_0)*sin(Omega_0));
y0 = r0*(cos(u)*sin(Omega_0) + sin(u)*cos(i_0)*cos(Omega_0));
z0 = r0*sin(u)*sin(i_0);
Vn = sqrt(mu/p)*(1+e*cos(thetta_0));
Vr = sqrt(mu/p)*e*sin(thetta_0);
vx0 = Vr*(cos(u)*cos(Omega_0) - sin(u)*cos(i_0)*sin(Omega_0)) - Vn*(cos(Omega_0)*sin(u) + sin(Omega_0)*cos(u)*cos(i_0));
vy0 = Vr*(cos(u)*sin(Omega_0) + sin(u)*cos(i_0)*cos(Omega_0)) - Vn*(sin(Omega_0)*sin(u) - cos(Omega_0)*cos(u)*cos(i_0));
vz0 = Vr*sin(u)*sin(i_0) + Vn*cos(u)*sin(i_0);
SC = [x0, y0, z0, vx0, vy0, vz0];
end
