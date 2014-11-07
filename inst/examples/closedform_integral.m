syms x y

mu = [0, 0];
sig = [1, 0; 0, 1];
%mvnpdf([0, 0], mu, sig)

lx = int(mvnpdf([x, y], mu, sig), x, -1, 1);
l = int(lx, y, -1, 1)