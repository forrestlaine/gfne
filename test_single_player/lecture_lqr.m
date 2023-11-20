%% Lecture LQR
close all; clear; clc;
T = 100;
n = 4;
m = 2;
dt = 0.05;
A = eye(n);
A(1,3) = dt;
A(2,4) = dt;
B = zeros(n,m);
B(1,1) = 0.5*dt*dt;
B(2,2) = 0.5*dt*dt;
B(3,1) = dt;
B(4,2) = dt;

c = zeros(n,1);
c(3) = 1;
c(4) = -0;

x0 = [-5;-5;0;2];
xf = [-3;3;0;0];

Q = 10*eye(n);
R = 100*eye(m);

q = zeros(n,1);
r = zeros(m,1);

blocks{1} = Q;
vecs{1} = q;
ind = 2;
for t = 1:T-1
    blocks{ind} = R;
    blocks{ind+1} = Q;
    vecs{ind} = r;
    vecs{ind+1} = q;
    ind = ind+2;
end
bigQ = blkdiag(blocks{:});
bigq = vertcat(vecs{:});

% bigA = zeros(T*n,T*n+(T-1)*m);
% bigb = zeros(T*n,1);
% bigA(1:n,1:n) = eye(n);
% bigb(1:n) = x0;
% ind = 0;
% for t = 1:T-1
%     bigA(t*n+1:t*n+n, ind+1:ind+n) = -A;
%     bigA(t*n+1:t*n+n, ind+n+1:ind+n+m) = -B;
%     bigA(t*n+1:t*n+n, ind+n+m+1:ind+2*n+m) = eye(n);
%     bigb(t*n+1:t*n+n) = c;
%     ind = ind+n+m;
% end

bigA = zeros(T*n+n,T*n+(T-1)*m);
bigb = zeros(T*n+n,1);
bigA(1:n,1:n) = eye(n);
bigb(1:n) = x0;
ind = 0;
for t = 1:T-1
    bigA(t*n+1:t*n+n, ind+1:ind+n) = -A;
    bigA(t*n+1:t*n+n, ind+n+1:ind+n+m) = -B;
    bigA(t*n+1:t*n+n, ind+n+m+1:ind+2*n+m) = eye(n);
    bigb(t*n+1:t*n+n) = c;
    ind = ind+n+m;
end
bigA(T*n+1:T*n+n,(T/2)*(n+m)+1:(T/2)*(n+m)+n) = eye(n);
bigb(T*n+1:T*n+n) = xf;

rk = rank(bigA);
[U,S,V] = svd(bigA);

V2 = V(:,rk+1:end);

pseudo_inv = bigA\bigb;
big_x_opt = pseudo_inv - V2/(V2'*bigQ*V2)*V2'*(bigQ*pseudo_inv+bigq);

xx = zeros(n,T);
uu = zeros(m,T-1);

ind = 0;
for t = 1:T-1
    xx(:,t) = big_x_opt(ind+1:ind+n);
    ind = ind+n;
    uu(:,t) = big_x_opt(ind+1:ind+m);
    ind = ind+m;
end
xx(:,T) = big_x_opt(ind+1:ind+n);
%%
figure; hold on;
for t = 1:T
    p = plot(xx(1,t), xx(2,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    axis([-10,10,-10,10]);
    pause(dt);
end