n = size(z0,1);
H = sparse(evalH(z0));
g = full(evalG(z0));
r = symrcm(H);

tic;
for i = 1:100
    z = randn(n,1);
    H = sparse(evalH(z));
    g = full(evalG(z));
    H(r,r)\g(r);
end
toc / 100