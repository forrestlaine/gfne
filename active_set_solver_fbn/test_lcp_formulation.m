H = randn(5,5);
H = H'*H;
f = randn(5,1);


A = randn(2,5);
b = randn(2,1);

[X,fval] = quadprog(H,f,-A,b);

M = [H -H; -H H];
q = [f; -f];
A = [A -A];

M = [M -A';
     A zeros(2,2)];
q = [q;b];

X2 = LCP(M,q);
X3 = X2(1:5)-X2(6:10);
D1 = [X2(1:5) X2(6:10)];

M2 = [H .1*eye(5)-H; .1*eye(5)-H H];
M2 = [M2 -A';
     A zeros(2,2)];

X4 = LCP(M2,q);
X5 = X4(1:5)-X4(6:10);
D2 = [X4(1:5) X4(6:10)];
