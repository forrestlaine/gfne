syms k1 k2 b3 b4
B = [0 0; 0 0; 0 0; b3 -b3; 0 0; b4 b4];
T = [.5 .5; -.5 .5];

Kh = [0 0 0 k1 0 0; 0 0 0 0 0 k2];

Bh = B*T;

Bh*Kh;

B*T*Kh

K = T*Kh

B*K