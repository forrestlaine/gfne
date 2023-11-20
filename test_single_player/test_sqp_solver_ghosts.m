%% test_sqp_solver

clear; clc; close all;
global linear_dyn;
linear_dyn = true;
global xga xgb xgc;
global half_border;
half_border = 10;
global sq_dist;

min_dist = 2;
sq_dist = min_dist*min_dist;

m{1} = 2;
m{2} = 2;
m{3} = 2;
N = 3;
n = 12;
T = 50;

x0a = [-3;-2;0;0];
x0b = [3;-3;0;0];
x0c = [0;3;0;0];
xga = [3;3;0;2];
xgb = [-3;3;0;2];
xgc = [0;-3;0;2];

x0 = [x0a;x0b;x0c];

l = cell(T+1,2);
h = cell(T+1,2);

for t = 1:T
    l{t,1} = @running_cost_a;
    l{t,2} = @running_cost_b;
    l{t,3} = @running_cost_c;
    
    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;
    
    
    g{t,1} = @inequality_constraint_a;
    g{t,2} = @inequality_constraint_b;
    g{t,3} = @inequality_constraint_c;

end
l{T+1,1} = @terminal_cost_a;
l{T+1,2} = @terminal_cost_b;
l{T+1,3} = @terminal_cost_c;
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;

h{T+1,1} = @terminal_constraint_a;
h{T+1,2} = @terminal_constraint_b;
h{T+1,3} = @terminal_constraint_c;
% h{T+1,1} = @empty_constraint;
% h{T+1,2} = @empty_constraint;

g{T+1,1} = @inequality_constraint_a;
g{T+1,2} = @inequality_constraint_b;
g{T+1,1} = @empty_constraint;
g{T+1,2} = @empty_constraint;
g{T+1,3} = @empty_constraint;

[residuals, X] = sqp_solver(@f, @avoid_constraint, h, g, l, n, m, N, T, x0);
%% plot solution
close all;
figure;
hold on;
axis([-2*half_border 2*half_border -2*half_border 2*half_border]);
drawrectangle('Position',[-5,-5,10,10]);
axis('equal');
for t = 1:T+1
    if linear_dyn
        ext_a = viscircles([X(1,t),X(2,t)],min_dist/2,'Color','b');
        ext_b = viscircles([X(5,t),X(6,t)],min_dist/2,'Color','r');
        ext_c = viscircles([X(9,t),X(10,t)],min_dist/2,'Color','g');
        path_a = plot(X(1,t),X(2,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
        path_b = plot(X(5,t),X(6,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'r');
        path_c = plot(X(9,t),X(10,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'g');
    else
        path_a = quiver(X(1,t),X(2,t),X(3,t)*cos(X(4,t)),X(3,t)*sin(X(4,t)),0,'b','linewidth',3,'MaxHeadSize',1);
        path_b = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
        path_c = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
    end
    axis([-10 10 -10 10]);
    pause(0.2);
    delete(path_a);
    delete(path_b);
    delete(path_c);
    delete(ext_a);
    delete(ext_b);
    delete(ext_c);
end



%%

function const = avoid_constraint(x)
    const = zeros(0,1);
%     xa = x(1:2);
%     xb = x(5:6);
%     xc = x(9:10);
%     const_ab = (xa-xb)'*(xa-xb)-30;
%     const_ac = (xa-xc)'*(xa-xc)-30;
%     const_bc = (xb-xc)'*(xb-xc)-30;
%     const = [const_ab;
%              const_ac;
%              const_bc];
end


function const = empty_constraint(~)
    const = zeros(0,1);
end


function val = empty_cost(~)
    val = 0;
end

function const = inequality_constraint_a(x)
    global half_border;
    global sq_dist;
    xa = x(1);
    ya = x(2);
    const = [xa + half_border;
             -xa + half_border;
             ya + half_border;
             -ya + half_border];
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
    const_ab = (xa-xb)'*(xa-xb)-sq_dist;
    const_ac = (xa-xc)'*(xa-xc)-sq_dist;
    const = [const;
             const_ac];
end

function const = inequality_constraint_b(x)
    global half_border;
    global sq_dist;
    xb = x(5);
    yb = x(6);
    const = [xb + half_border;
            -xb + half_border;
             yb + half_border;
            -yb + half_border];
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
    const_ab = (xa-xb)'*(xa-xb)-sq_dist;
    const_bc = (xb-xc)'*(xb-xc)-sq_dist;
    const = [const;
             const_ab;
             const_bc];
end

function const = inequality_constraint_c(x)
    global half_border;
    global sq_dist;
    xc = x(9);
    yc = x(10);
    const = [xc + half_border;
            -xc + half_border;
             yc + half_border;
            -yc + half_border];
end

function const = terminal_constraint_a(x)
    const = zeros(0,1);
    global xga;
    xa = x(1:4);
    const = xa-xga;
end

function const = terminal_constraint_b(x)
    const = zeros(0,1);
    global xgb;
    xb = x(5:8);
    const = xb-xgb;
end

function const = terminal_constraint_c(x)
    const = zeros(0,1);
    global xgc;
    xc = x(9:12);
    const = xc-xgc;
end


function val = terminal_cost_a(x)
    val = 0;
end

function val = terminal_cost_b(x)
    val = 0;
end

function val = terminal_cost_c(x)
    val = 0;
end

function val = running_cost_a(xu)
    ua = xu(13:14);
    val = ua'*ua;
end

function val = running_cost_b(xu)
    ub = xu(15:16);
    val = ub'*ub;
end

function val = running_cost_c(xu)
    uc = xu(17:18);
    val = uc'*uc;
end

function vec = dubins_dyn(x,u)
    global linear_dyn;

    if linear_dyn
        vec = x(3);
        vec = [vec; x(4)];
        vec = [vec; u(1)];
        vec = [vec; u(2)];
    else
        vec = x(3)*cos(x(4));
        vec = [vec; x(3)*sin(x(4))];
        vec = [vec; u(1)];
        vec = [vec; u(2)];
    end
end

function next = f(xu)
    xa = xu(1:4);
    xb = xu(5:8);
    xc = xu(9:12);
    ua = xu(13:14);
    ub = xu(15:16);
    uc = xu(17:18);
    
    za = xa + 0.1*dubins_dyn(xa,ua);
    zb = xb + 0.1*dubins_dyn(xb,ub);
    zc = xc + 0.1*dubins_dyn(xc,uc);
    % maybe wrap angles( if it doesn't break differentiability)
    next = [za;zb;zc];
end