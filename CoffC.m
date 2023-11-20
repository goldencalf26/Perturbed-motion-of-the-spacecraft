function [Cxx, Cyy, Czz] = CoffC(X)
Cxx = X(:,2).*X(:,6) - X(:,3).*X(:,5);
Cyy = X(:,3).*X(:,4) - X(:,1).*X(:,6);
Czz = X(:,1).*X(:,5) - X(:,2).*X(:,4);
end
