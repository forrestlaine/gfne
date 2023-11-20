clear; clc; close all;
import casadi.*

x = SX.sym('x'); 
px = SX.sym('px');
py = SX.sym('py');
% y = SX.sym('y');

s = 0.2;
h = 3;
w = 1.5;

pxc = 2;
pyc = 3;

y1 = -s*x + h/w^2 * (w^4-x.^4).^(1/4);
y2 = -s*x - h/w^2 * (w^4-x.^4).^(1/4);
d1 = (x-pxc)^2 + (y1-pyc)^2;
d2 = (x-pxc)^2 + (y2-pyc)^2;

gd1 = gradient(d1,x);
gd2 = gradient(d2,x);

evaly1 = Function('y1',{x},{y1});
evaly2 = Function('y2',{x},{y2});
evald1 = Function('d1',{x,px,py},{d1});
evald2 = Function('d2',{x,px,py},{d2});
evalg1 = Function('gd1', {x},{gd1});
evalg2 = Function('gd2', {x},{gd2});

G = rootfinder('G','newton',evalg1);

N = 50;
XX = linspace(-w,w,N);


YY1 = full(evaly1(XX));
YY2 = full(evaly2(XX));
D1 = full(evald1(XX,pxc,pyc));
D2 = full(evald2(XX,pxc,pyc));
GD1 = full(evalg1(XX));
GD2 = full(evalg2(XX));

figure;
hold on;
plot(XX,YY1);
plot(XX,YY2);
plot(XX,GD1);
plot(XX,GD2);
axis('equal')





% x = linspace(-1,1, 100);
% xx = linspace(-0.99,0.99,100);
% y1 = -s*x + h/w^2 * (w^4-x.^4).^(1/4);
% dy1 = -s - h/w^2 * xx.^3 ./ (w^4 -xx.^4).^(3/4);
% 
% y2 = -s*x - h/w^2 * (w^4-x.^4).^(1/4);
% dy2 = -s + h/w^2 * xx.^3 ./ (w^4 -xx.^4).^(3/4);
% figure(); hold on;
% plot(x,y1);
% plot(xx,dy1);
% plot(x,y2);
% plot(xx,dy2);
% axis('equal');
% 
% % plot(x,y2);
% % axis([-5 5 -5 5])