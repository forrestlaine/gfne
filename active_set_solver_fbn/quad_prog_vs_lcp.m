n = 5;
p = 3;


Q = randn(n,n);
Q = Q'*Q;
% Q = zeros(n,n);

q = randn(n,1);


A = randn(p,n);
b = randn(p,1);%-10*ones(p,1);
% b = [randn(5,1)+3*ones(n,1); randn(5,1)+3*ones(n,1)];
     
% p = 10;
%%
L = -100*ones(n,1);


X = quadprog(Q,q,-A,b);

M = [Q -A'; A zeros(p)];
q = [q; b];
l = [-inf*ones(n,1); zeros(p,1)];
u = inf*ones(n+p,1);

[z,retcode] = LMCP(M,q,l,u);


x = z(1:n);

x-X
