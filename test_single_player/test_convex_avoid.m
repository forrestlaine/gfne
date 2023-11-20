clear; clc; close all;
N = 100;

global X; 
global Y;
global half_w;
global half_l;
global shear;
half_w = 1.5;
half_l = 5.0;
shear = 0.2;
X = [1.3;-1.2];
Y = [0;1.0];

S = fmincon(@shear_objective, [0,0,0,0], [],[],[],[],[],[],@shear_constraint);
shear_constraint([S(3),S(4),S(1),S(2)])
figure; hold on;
plot(S(1),S(2),'or','MarkerSize', 6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'r');
plot(S(3),S(4),'or','MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');
[lat,long] = drawShearedRecircle(X(2),X(1),shear,half_w,half_l);
plot(lat,long,'-b');
[lat,long] = drawShearedRecircle(Y(2),Y(1),shear,half_w,half_l);
plot(lat,long,'-r');
plot([-1,1],[0,0],'-k');
plot([0,0],[-1,1],'-k');

axis('equal')
function val = shear_objective(x)
    global X Y half_w half_l shear;
    pa = [x(1);x(2)];
    pb = [x(3);x(4)];
    val = 0;
    val = 0.1*(pa-pb)'*(pa-pb);
    val = val + ((pa(1)-Y(1))/half_w)^4 + (((pa(2)-Y(2)) + shear*(pa(1)-Y(1)))/half_l)^4;
    val = val + ((pb(1)-X(1))/half_w)^4 + (((pb(2)-X(2)) + shear*(pb(1)-X(1)))/half_l)^4;
    val = val*10000;
end

function [C,Ceq] = shear_constraint(x)
    half_w = 1.5;
    half_l = 5.0;
    shear = 0.2;
    global X Y;
    pa = [x(1);x(2)];
    pb = [x(3);x(4)];
    C = [((pa(1)-X(1))/half_w)^4 + (((pa(2)-X(2)) + shear*(pa(1)-X(1)))/half_l)^4 - 1;
         ((pb(1)-Y(1))/half_w)^4 + (((pb(2)-Y(2)) + shear*(pb(1)-Y(1)))/half_l)^4 - 1];
    Ceq = [];
end





