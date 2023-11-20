% syms k l g m
% A = [1 2 1;
%     2 4 2;
%     4 8 0;
%     8 1 1];
% y = [2; 4; 5];
% rank(A)
% 
% x = linspace(-3,3,100);
% y1 = sin(x)+1;
% y2 = 2*exp(x);
% 
% figure;
% hold on;
% plot(x,y1);
% plot(x,y2);

A = [0 1; -1 -2.5];
[V D] = eig(A)
T = 0.2;
V = [-2 1; 1 -2];


Ad = V*[exp(-0.5*T) 0; 0 exp(-2*T)]*inv(V)
Bd = V*[2/3*(exp(-T/2) - 1); 1/3*(exp(-2*T) -1)]