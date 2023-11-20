clear; clc; close all;

H = randn(3,3);
H = H'*H + eye(3);
h = randn(3,1);


A = randn(1,3);
B = randn(1,3);

a = randn(1);
b = randn(1);

sk = 0.1;
zk = 0.33;
yk = randn(1);
mu = 0.1;

Sig = zk/sk;

MM = [H zeros(3,1) A' B';
      zeros(1,3), Sig, zeros(1,1) -1;
      A 0 0 0;
      B -1 0 0];
  
mm = [-h + A'*yk + B'*zk;
     mu/sk-zk;
     -a;
     sk-b];
 
sol = MM\mm;
px = sol(1:3);
ps = sol(4);
py = -sol(5);
pz = -sol(6);

M2 = MM;
m2 = [-h; mu/sk; -a; sk-b];

sol2 = MM\m2;
px2 = sol2(1:3);
ps2 = sol2(4);
py2 = -sol2(5)-yk;
pz2 = -sol2(6)-zk;

M3 = [H A' B'; A zeros(1,2); B 0 -1/Sig];
m3 = [-h; -a; sk-b+mu/zk];

sol3 = M3\m3;
px3 = sol3(1:3)
py3 = -sol3(4)-yk;
pz3 = -sol3(5)-zk;

M4 = [H+B'*Sig*B A'; A 0];
m4 = [-h - B'*Sig*b + B'*zk + B'*mu/sk; 
      -a];

sol4 = M4\m4;
px4 = sol4(1:3)
