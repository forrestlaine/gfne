clear;
clc;
close all;

Q1 = [2 -1; -1 1];
q1 = [-1;0];
zer = -Q1\q1;


A = [-.85, 1];
mult = 2;
A = mult*A;

M = [Q1, A'; A, 0];

h0 = -1*mult;
h1 = 0* mult;

m = [-q1; h0];
m2 = [-q1; h1];

soln = M\m;
opt = soln(1:2);
lam = soln(3);
soln2 = M\m2;
opt2 = soln2(1:2);
lam2 = soln2(3);

zerb = [-.5; .5];
Q2 = eye(2);
q2 = -zerb;

v0 = 0.5*opt2'*Q1*opt2 + opt2'*q1;
v1 = 0.5*opt'*Q1*opt + opt'*q1;

% dcost_a = (v1-v0) / h1

% lam

v0b = 0.5*opt2'*Q2*opt2 + opt2'*q2;
v1b = 0.5*opt'*Q2*opt + opt'*q2;

dcost_b = (v1b-v0b) / (h1-h0);

dc = (v0b - v1b);

gx2 = Q2*opt + q2;
grad1 = Q1*opt + q1;
gx1 = (A/Q1)';


X = linspace(-3,3,500);
Y = linspace(-2,6,500);
for i = 1:500
    for j = 1:500
        vec = [X(i);Y(j)];
        Z1(i,j) = 0.5*vec'*Q1*vec + vec'*q1;
        Z2(i,j) = 0.5*vec'*Q2*vec + vec'*q2;
        
    end
end
figure;
hold on;

contour(X,Y,Z1',[v0,v1],'r', 'DisplayName', 'Cost Contour 1');
contour(X,Y,Z2',[v0b, v1b],'b', 'DisplayName', 'Cost Contour 2');

% Quadratic zeros
plot(zer(1), zer(2),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
plot(zerb(1), zerb(2),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

% Optimal spot for h0
plot(opt(1), opt(2),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
% Optimal spot for h0
plot(opt2(1), opt2(2),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'g');



bnd = 2;

plot([-bnd,bnd],[(A(1)*bnd+h0)/A(2),(-A(1)*bnd+h0)/A(2)], 'g')
plot([-bnd,bnd],[(A(1)*bnd+h1)/A(2),(-A(1)*bnd+h1)/A(2)], 'g')

dd = (A/Q1*A')\A/Q1;
slam = -dd*gx2

p1 = plot([opt(1), opt(1)-gx2(1)], [opt(2), opt(2)-gx2(2)],'--b','linewidth',1, 'DisplayName', 'Cost Gradient 1');
p2 = plot([opt(1), opt(1)-grad1(1)], [opt(2), opt(2)-grad1(2)],'--r','linewidth',1, 'DisplayName', 'Cost Gradient 2');
p3 = plot([opt(1), opt(1)-A(1)], [opt(2), opt(2)-A(2)],'--g','linewidth',1, 'DisplayName', 'Constraint Gradient');
% plot([opt(1), opt(1)+A(1)], [opt(2), opt(2)+A(2)],'--k','linewidth',1, 'DisplayName', 'Cost Gradient 1');
p4 = plot([opt(1), opt(1)+gx1(1)], [opt(2), opt(2)+gx1(2)],'--k','linewidth',1, 'DisplayName', 'Constrained Path');
plot([opt(1), opt(1)-gx1(2)], [opt(2), opt(2)+gx1(1)],'-.k','linewidth',1);
plot([opt(1), opt(1)+gx1(2)], [opt(2), opt(2)-gx1(1)],'-.k','linewidth',1);
plot([opt(1), opt(1)+dd(1)], [opt(2), opt(2)+dd(2)],':k','linewidth',1);
plot([opt(1)-gx2(1) - slam*A(1), opt(1)-gx2(1)], [opt(2)-gx2(2)-slam*A(2), opt(2)-gx2(2)],':k','linewidth',1);
% 
legend([p1,p2,p3, p4])

axis('equal')
axis([-2 2 -2 2])
% axis([.5-10*h1 .5+10*h1 -10*h1 10*h1])
