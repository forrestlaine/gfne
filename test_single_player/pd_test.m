n = 3;
m = 2;
l = 2;

Q = zeros(n);
B = rand(n,m);
G = rand(l,m);
I = eye(n);
Z = zeros(l,n);
R = eye(m);

M = [Q Z' I;
    Z -G/R*G' -G/R*B';
    I -B/R*G' -B/R*B'];

eig(M)