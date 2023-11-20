%% test_sqp_solver

clear; clc; close all;
global linear_dyn;
global xga;
global xgb;
linear_dyn = true;

m{1} = 2;
m{2} = 2;
N = 2;
n = 8;
T = 50;

if linear_dyn
    x0a = [-3;-2;0;0];
    x0b = [3;-3;0;0];
    xga = [3;3;0;5];
    xgb = [-3;3;0;5];
else
    x0a = [-3;-2;0;pi/2];
    x0b = [3;-3;0;pi/2];
    xga = [3;3;1;pi/2];
    xgb = [-3;3;1;pi/2];
end
x0 = [x0a;x0b];

l = cell(T+1,2);
h = cell(T+1,2);

for t = 1:T
    l{t,1} = @running_cost_a;
    l{t,2} = @running_cost_b;
    
    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    
    
    g{t,1} = @inequality_constraint_a;
    g{t,2} = @inequality_constraint_b;
%     g{t,1} = @empty_constraint; 
%     g{t,2} = @empty_constraint; 
    

end
l{T+1,1} = @terminal_cost_a;
l{T+1,2} = @terminal_cost_b;
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
h{T+1,1} = @terminal_constraint_a;
h{T+1,2} = @terminal_constraint_b;
% h{T+1,1} = @empty_constraint;
% h{T+1,2} = @empty_constraint;

g{T+1,1} = @inequality_constraint_a;
g{T+1,2} = @inequality_constraint_b;
g{T+1,1} = @empty_constraint;
g{T+1,2} = @empty_constraint;

rr{1} = false;
rr{2} = false;



[residuals, X] = sqp_solver(@f, @avoid_constraint, h, g, l, n, m, N, T, rr, x0);
%% plot solution
close all;
figure;
hold on;
axis([-10 10 -10 10]);
drawrectangle('Position',[-5,-5,10,10]);
axis('equal');
for t = 1:T+1
    if linear_dyn
        path_a = plot(X(1,t),X(2,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        path_b = plot(X(5,t),X(6,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    else
        path_a = quiver(X(1,t),X(2,t),X(3,t)*cos(X(4,t)),X(3,t)*sin(X(4,t)),0,'b','linewidth',3,'MaxHeadSize',1);
        path_b = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
    end
    axis([-10 10 -10 10]);
    pause(0.1);
    delete(path_a);
    delete(path_b);
end



%%

function const = avoid_constraint(x)
    const = zeros(0,1);
    xa = x(1:2);
    xb = x(5:6);
    const = (xa-xb)'*(xa-xb)-30;
end


function const = empty_constraint(~)
    const = zeros(0,1);
end


function val = empty_cost(~)
    val = 0;
end

function const = inequality_constraint_a(x)
    xa = x(1);
    ya = x(2);
%     va = x(3);
    const = [xa + 5;
             -xa + 5;
             ya + 5;
             -ya + 5];
end

function const = inequality_constraint_b(x)
    xb = x(5);
    yb = x(6);
%     vb = x(7);
    const = [xb + 5;
            -xb + 5;
             yb + 5;
            -yb + 5];
end



function val = midpoint_cost_a(x)
    goal = [5;0];
    xa = x(1:2);
    const = xa-goal;
    val = 10*const'*const + running_cost_a(x);
end

function val = midpoint_cost_b(x)
    goal = [-5;0];
    xb = x(5:6);
    const = xb-goal;
    val = 10*const'*const + running_cost_b(x);
end

function const = terminal_constraint_a(x)
    global xga;
    xa = x(1:4);
    const = xa-xga;
end

function const = terminal_constraint_b(x)
    global xgb;
    xb = x(5:8);
    const = xb-xgb;
end

function val = terminal_cost(~)
    val = 0;
end

function val = terminal_cost_a(x)
    val = 0;
    goal = [3;3;0;1];
    xa = x(1:4);
    xb = x(5:8);
    val = 100*(xa-goal)'*(xa-goal);
end

function val = terminal_cost_b(x)
    val = 0;
    goal = [-3;3;0;1];
    xa = x(1:4);
    xb = x(5:8);
    val = 100*(xb-goal)'*(xb-goal);
end

function val = running_cost_a(xu)
    ua = xu(end-4+1:end-4+2);
    val = 1*ua'*ua;
    xa = xu(1:4);
    val = val;
end

function val = running_cost_b(xu)
    ub = xu(end-2+1:end);
    val = 1*ub'*ub;
    xb = xu(5:8);
    val = val;
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
    ua = xu(9:10);
    ub = xu(11:12);
    
    za = xa + 0.1*dubins_dyn(xa,ua);
    zb = xb + 0.1*dubins_dyn(xb,ub);
    % maybe wrap angles( if it doesn't break differentiability)
    next = [za;zb];
end