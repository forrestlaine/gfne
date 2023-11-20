clear; clc;

A = eye(8); 
A(1,3) = 0.1;
A(2,4) = 0.1;
A(5,7) = 0.1;
A(6,8) = 0.1;

B1 = [0 0; 0 0; 0.1 0; 0 0.1];
B2 = B1;

B = [B1 zeros(4,2); zeros(4,2) B2];


Q1 = blkdiag(zeros(4),eye(4));
Q2 = [eye(4) -eye(4); -eye(4) eye(4)];
% Q1 = zeros(8);
% Q2 = zeros(8);

G1 = zeros(0,8);
G2 = zeros(0,8);
% G1 = [zeros(4) eye(4)];
% G2 = [-eye(4) eye(4)];

x0 = [0 0 0 0 1 1 -10 -10]';

M2 = [eye(4) zeros(4,8) B' zeros(4,16);
     zeros(8,4) Q1 -eye(8) zeros(8) G1' zeros(8,4);
     zeros(8,4) Q2 zeros(8) -eye(8) zeros(8,4) G2';
     B -eye(8) zeros(8,24);
     zeros(8,4) [G1;G2] zeros(8,24)];

Mx = [zeros(4,8); zeros(8); zeros(8); -A; zeros(8)];

K = pinv(M2)*Mx;

sol = K*x0;
u = sol(1:4);

G2*(A*x0 + B*u)
G1*(A*x0+B*u)

Ma = [eye(2) zeros(2,8) B1' zeros(2,4) zeros(2,4);
      zeros(8,2) Q1 -eye(8) G1';
      [B1; zeros(4,2)] -eye(8) zeros(8,12);
      zeros(4,2) G1 zeros(4,12)];
  
Nau = [zeros(2); zeros(8,2); zeros(4,2); B2; zeros(4,2)];
Nax = [zeros(2,8); zeros(8); -A; zeros(4,8)];

Mb = [eye(2) zeros(2,8) zeros(2,4) B2' zeros(2,4);
      zeros(8,2) Q2 -eye(8) G2';
      [zeros(4,2); B2] -eye(8) zeros(8,12);
      zeros(4,2) G2 zeros(4,12)];
  
Nbu = [zeros(2); zeros(8,2); B1; zeros(4,2); zeros(4,2)];
Nbx = [zeros(2,8); zeros(8); -A; zeros(4,8)];

Kau = pinv(Ma)*Nau;
Kax = pinv(Ma)*Nax;

Kau = Kau(1:2,:);
Kax = Kax(1:2,:);

Kbu = pinv(Mb)*Nbu;
Kbx = pinv(Mb)*Nbx;

Kbu = Kbu(1:2,:);
Kbx = Kbx(1:2,:);

MM = [eye(2) Kau; Kbu eye(2)];
MMx = [Kax; Kbx];
K2 = MM\MMx;

K1 = K(1:4,:);
F1 = A+B*K1;
F2 = A+B*K2;

Z1 = [eye(8); K1];
Z2 = [eye(8); K2];

P1a = Z1'*blkdiag(Q1,eye(2),zeros(2))*Z1 + F1'*Q1*F1;
P1b = Z2'*blkdiag(Q1,eye(2),zeros(2))*Z2 + F2'*Q1*F2;

P2a = Z1'*blkdiag(Q2,zeros(2),eye(2))*Z1 + F1'*Q1*F1;
P2b = Z2'*blkdiag(Q2,zeros(2),eye(2))*Z2 + F2'*Q2*F2;


x0'*(P1a-P1b)*x0
x1 = [0 0 0 0 -.5 -.5 5 5]';
K1*x1

x0'*(P2a-P2b)*x0
x1'*(P2a-P2b)*x1
 
P2a-P2b

% (A*x0 + B*[10;10;-10;-10])
