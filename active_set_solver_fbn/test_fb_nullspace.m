% test_fb_nullspace

A = randn(8,8);
g = randn(2,1);
c = [g; zeros(6,1)];

AA = [A c; c' 0];

B = randn(8,3);
BB = [B; zeros(1,3)];

K = A\B;
Ka = K(1:2,:);
Kb = K(3:4,:);

Ai = inv(A);
Ail = Ai(:,1:2);
ggt = g*g';
M = Ail*ggt / (c'*Ai*c);
Mu = M(3:4,:);

KK = Kb - Mu*Ka;

KK2 = AA\BB;
KK2 = KK2(3:4,:);