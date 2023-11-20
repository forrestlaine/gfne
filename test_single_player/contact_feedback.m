%%
clear; clc; close all;
%% Problem setup

Q = 1*eye(2);
R = 1*eye(1);

T = 0.1;
B = [T*T/2; T];
C = B;
A = [1 T; 0 1];

L = [2/(T*T) 2/T];

%% 

K1 = (R+B'*Q*B)\(B'*Q*A);

Ah = (A-C*L);
Bh = (B-C);

K2 = (R+Bh'*Q*Bh)\(Bh'*Q*Ah);
%%

MM = [Q+A'*Q*A A'*Q*C A'*Q*B;
       C'*Q*A C'*Q*C C'*Q*B;
       B'*Q*A B'*Q*C R+B'*Q*B];
Kx = (R+B'*Q*B)\(B'*Q*A);
Kl = (R+B'*Q*B)\(B'*Q*C);

D = [eye(2) [0;0];
     0 0 1;
     -Kx -Kl];
 
 V = D'*MM*D;
 
 M3 = A-B*Kx;
 M3 = M3(1,:);
 
 V1 = V(1:2,1:2);
 
 S1 = A-B*Kx;
 S2 = C-B*Kl;
 S1 = S1(1,:);
 S2 = S2(1);
 G = S1/S2;
 
 D2 = [eye(2);
       -G];
 V2 = D2'*V*D2;
 
 % cond3*X > 0;



%%

% Condition 1: (L-K1)x >= 0
M1 = L-K1
% Condition 2: (L-K2)x < 0
M2 = L-K2

M3

bx = 100;
by = 100;


X = linspace(-bx,bx,100);
Y = linspace(-by,by,100);
Y1 = -M1(1)/M1(2)*X;
Y2 = -M2(1)/M2(2)*X;
Y3 = -M3(1)/M3(2)*X;

Z = zeros(100,100);
for i = 1:100
    for j = 1:100
        s = [X(i);Y(j)];
        if (M3*s > 0)
            Z(i,j) = 100*0.5*s'*V1*s;
        else
            Z(i,j) = 0.5*s'*V2*s;
        end
    end
end


%%
figure; hold on;
% plot(X,Y1, 'linewidth', 3);
% plot(X,Y2, 'linewidth', 3);
plot(X,Y3, 'linewidth', 3);

% figure; hold on;
contour(X,Y,Z',100)
axis([-bx bx -by by])

% x_bounds = [-bx bx bx 0 -bx];
% y1_bounds = [max(-M1(1)/M1(2)*-bx, by), max(-M2(1)/M2(2)*bx, by), max(-M2(1)/M2(2)*bx, -by), 0, max(-M1(1)/M1(2)*-bx, -by)];
% y2_bounds = [min(-M2(1)/M2(2)*-bx, by), 0, min(-M1(1)/M1(2)*bx, by), min(-M1(1)/M1(2)*bx, -by), min(-M2(1)/M2(2)*-bx, -by)];
% 
% fill(x_bounds, y1_bounds, 'b');
% fill(x_bounds, fliplr(y2_bounds), 'r');
% fill([0, bx, bx], [0, -M1(1)/M1(2)*bx,-M2(1)/M2(2)*bx], 'm');
