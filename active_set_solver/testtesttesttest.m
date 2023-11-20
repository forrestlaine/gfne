% test px py K matrices

R11 = eye(2);
R12 = randn(2,2);
R21 = randn(2,2);
R22 = eye(2);

Q1 = eye(3);
Q2 = eye(3);

B1 = randn(3,2);
B2 = randn(3,2);
G1 = randn(1,3);
G2 = randn(1,3);

M1 = blkdiag([R11 R12; R21 R22], [Q1;Q2]);
M2 = [B1' zeros(2,3); zeros(2,3) B2'; -eye(3) zeros(3); zeros(3) -eye(3)];
M3 = [zeros(4,2); G1' zeros(3,1); zeros(3,1) G2'];
M4 = [B1 B2 -eye(3); zeros(1,4) G1; zeros(1,4) G2];
M5 = zeros(size(M4,1),size([M2 M3],2));
M = [M1 M2 M3; M4 M5];
Mi = inv(M);


