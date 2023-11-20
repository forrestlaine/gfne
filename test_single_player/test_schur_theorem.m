n = 6;
m = 3;

Q = randn(n,n);
Q = Q'*Q;
% Q = zeros(n,n);

S = randn(m,n);
S = zeros(m,n);

R = randn(m,m);
R = R'*R+0.1*eye(m);
R = 0.1*eye(m);

L = randn(1,n+m);
% L = zeros(1,n+m);
% L(5) = 1;
N = randn(1,n+m);
% N = zeros(1,n+m);
% N(5) = 1;
D = 2*eye(1);

Md = -N'*L-L'*N+L'*D*L;
eig(Md(n+1:end,n+1:end))
M1 = [Q S'; S R]-N'*L-L'*N+L'*D*L;
M2 = [Q S'; S R];

Q1 = M1(1:n,1:n);
S1 = M1(n+1:n+m,1:n);
R1 = M1(n+1:n+m,n+1:n+m);

Q2 = Q;
S2 = S;
R2 = R;

L1 = Q1-S1'/R1*S1;
L2 = Q2-S2'/R2*S2;

eig(M1-M2)
eig(L1-L2)

Lq = Md(1:n,1:n);
Ls = Md(n+1:n+m,1:n);
Lr = Md(n+1:n+m,n+1:n+m);

% diff = Lq - (S+Ls)'*(inv(R) - inv(R+R*inv(Lr)*R))*(S+Ls) + S'*inv(R)*S;
% diff2 = Lq + (S+Ls)'*inv(R+R*inv(Lr)*R)*(S+Ls) - Ls'*inv(R)*S -S'*inv(R)*Ls - Ls'*inv(R)*Ls;

% A = S+Ls;
% B = inv(R)*Ls;
% C = inv(R+R*inv(Lr)*R);
% eig(Lq+A'*C*A-B'*A-A'*B)


