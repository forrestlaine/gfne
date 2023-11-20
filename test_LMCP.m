n = 5;
p = 3;


Q = randn(4,n);
Q = Q'*Q;
% Q = zeros(n,n);

q = randn(n,1);


A = randn(p,n);
b = randn(p,1)-100*ones(p,1);
% b = [randn(5,1)+3*ones(n,1); randn(5,1)+3*ones(n,1)];
     
% p = 10;
%%
L = -1000*ones(n,1);
% L = -inf*ones(n,1);

[X,FVAL,EXITFLAG] = quadprog(Q,q,-A,b);
% if EXITFLAG ~= 1
%     disp('Infeasible!');
% end

M = [Q -A'; A zeros(p)];
q = [q; b];
l = [L; zeros(p,1)];
u = inf*ones(n+p,1);

[z,retcode] = LMCP(M,q,l,u);
[zz,rr] = LCPSolve(M,q+M*l);

x = z(1:n);

solution_error = norm(x-X)