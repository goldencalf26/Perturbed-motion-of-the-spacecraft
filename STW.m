function [S_w, T_w, W_w] = STW(eps, r_mod, inc_w, u_w);
S_w = eps./r_mod.^4.*(3.*sin(inc_w).^2.*sin(u_w).^2 - 1);
T_w = -eps./r_mod.^4.*sin(inc_w).^2.*sin(2.*u_w);
W_w = -eps./r_mod.^4.*sin(2.*inc_w).*sin(u_w);
end
