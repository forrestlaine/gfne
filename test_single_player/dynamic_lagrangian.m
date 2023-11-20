clear; clc; close all;

Q = diag([0,1,0]);
R = rand(2,2);
R = R'*R + eye(2);
QQ = blkdiag(Q,Q,Q,Q);
RR = blkdiag(R,R,R);
A = rand(3,3);
B = rand(3,2);
AA = [eye(3), zeros(3,3*3+3*2);
    -A eye(3) zeros(3,6) -B zeros(3,4);
    zeros(3,3), -A eye(3) zeros(3,5) -B, zeros(3,2);
    zeros(3,6), -A eye(3) zeros(3,4) -B];

m = 6;
n = 12;
BB = AA(:,n+1:end);
AA = AA(:,1:n);





MM = [RR, zeros(m,n), -BB';
    zeros(n,m), QQ, -AA';
    BB, AA, zeros(n,n)];

q = rand(3,1);
r = rand(2,1);
c = rand(3,1);
qq = [q;q;q;q];
rr = [r;r;r];
cc = [zeros(3,1); c; c; c];


AAi = inv(AA);
RRi = inv(RR);

nn = size(AAi,1);

MMi = inv(MM);
soln = MMi*[-rr;-qq;cc];
u = soln(1:m);
x = soln(m+1:m+n);
lams = soln(m+n+1:end);

e1 =             RRi*BB'/(AA'+QQ*AAi*BB*RRi*BB')*QQ*AAi;
e2 = AAi - AAi*BB*  RRi*BB'/(AA'+QQ*AAi*BB*RRi*BB')*QQ*AAi;
ee = (AA' + QQ*AAi*BB*RRi*BB')\[QQ*AAi*BB*RRi, eye(nn)];

dux_dh = [e1;e2];

lams = dux_dh'*([RR zeros(m,n); zeros(n,m) QQ]*[u;x] + [rr;qq]);
