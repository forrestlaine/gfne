
% B = [1; -1; 1];
% A = [1 0 0.1; 0 1 0.1; 0 0 1];
% 
% rank([B A*B A*A*B])
% 
% B = [-1; 1; 2];
% A = [1 0 0.1; 0 1 0.1; 0 0 1];
% 
% rank([B A*B A*A*B])
% 
% B = [2; 0; 3];
% A = [1 0 0.5; 0 1 0.75; 0 0 1];
% 
% rank([B A*B A*A*B])
% 
% eig([-2 -1; 3.75 2])
% eig([3 1; -8.75 -3])
% eig([2 1; -3.75 -2])

X = linspace(-2,4*pi, 100);
Y2 = 2*exp(-X);
Y1 = sin(X)+1;
figure;
hold on;
plot(X,Y1);
plot(X,Y2);
xline(0);
yline(0);
legend('f1(x1,x2) = 0', 'f2(x1,x2) = 0');
xlabel('X1');
ylabel('X2');