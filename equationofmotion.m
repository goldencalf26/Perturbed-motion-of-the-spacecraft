function dXdt = equationofmotion(t, X, sigma, ro, mu, J)
r = norm(X(1:3));
V = norm(X(4:6));
dx1 = X(4);
dx2 = X(5);
dx3 = X(6);
dx4 = -mu./r.^3.*X(1) + (3/2*J*mu./r.^2.*(6371./r).^2*(5.*X(3).^2/r.^2-1)).*X(1)./r - sigma*ro.*V.*X(4);
dx5 = -mu./r.^3.*X(2) + (3/2*J*mu./r.^2.*(6371./r).^2*(5.*X(3).^2/r.^2-1)).*X(2)./r - sigma*ro.*V.*X(5);
dx6 = -mu./r.^3.*X(3) + (3/2*J*mu./r.^2.*(6371./r).^2*(5.*X(3).^2/r.^2-3)).*X(3)./r - sigma*ro.*V.*X(6);
dXdt = [dx1;dx2;dx3;dx4;dx5;dx6];
%end
end



