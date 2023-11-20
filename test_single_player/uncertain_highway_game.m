%% test_sqp_solver

clear; clc; close all;
global linear_dyn;
linear_dyn = true;

global xga xgb xgc xgd;
global probc;
global min_dist sq_dist;
global border_left border_right;

border_left = -5;
border_right = 5;

probc = 0;

min_dist = 3;
car_rad = 1;

sq_dist = min_dist*min_dist;

m{1} = 2;
m{2} = 2;
m{3} = 2;
m{4} = 2;


N = 4;
n = 4*N;
T = 50;


x0a = [2.5;10;0;0];
x0b = [2.5;15;0;0];
x0c = [-2.5;-10;0;0];
x0d = [-2.5;-10;0;0];


xga = [-2.5;20;0;0];
xgb = [2.5;15;0;0];
xgc = [-2.5;10;0;0];
xgd = [-2.5;25;0;0];

x0 = [x0a;x0b;x0c;x0d];

l = cell(T+1,2);
h = cell(T+1,2);

for t = 1:T
    l{t,1} = @running_cost_a;
    l{t,2} = @running_cost_b;
    l{t,3} = @running_cost_c;
    l{t,4} = @running_cost_d;
    
    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;
    h{t,4} = @empty_constraint;
   
    g{t,1} = @inequality_constraint_a;
    g{t,2} = @inequality_constraint_b;
    g{t,3} = @inequality_constraint_c;
    g{t,4} = @inequality_constraint_d;
    
%     g{t,1} = @empty_constraint;
%     g{t,2} = @empty_constraint;
%     g{t,3} = @empty_constraint;
%     g{t,4} = @empty_constraint;

end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;
l{T+1,4} = @empty_cost;

h{T+1,1} = @terminal_constraint_a;
h{T+1,2} = @terminal_constraint_b;
h{T+1,3} = @terminal_constraint_c;
h{T+1,4} = @terminal_constraint_d;
% h{T+1,1} = @empty_constraint;
% h{T+1,2} = @empty_constraint;
h{T+1,3} = @empty_constraint;
h{T+1,4} = @empty_constraint;

% g{T+1,1} = @terminal_ineq_constraint_a;
% g{T+1,2} = @terminal_ineq_constraint_b;
% g{T+1,2} = @empty_constraint;
% g{T+1,3} = @terminal_ineq_constraint_c;
% g{T+1,4} = @terminal_ineq_constraint_d;

g{T+1,1} = @inequality_constraint_a;
g{T+1,2} = @inequality_constraint_b;
g{T+1,3} = @inequality_constraint_c;
g{T+1,4} = @inequality_constraint_d;

% g{T+1,3} = @empty_constraint;
% g{T+1,4} = @empty_constraint;


[X,z,duration,eps_opt] = ai_solver(@f, h, g, l, n, m, N, T, 1, x0);

for t=1:T
    l{t,1} = @running_cost_a_polite;
end

[Xp,z,durationpp,eps_optp] = ai_solver(@f, h, g, l, n, m, N, T, 1, x0,z);


%% plot solution
close all;
figure;
hold on;
axis([-10 10 -10 10]);
plot([border_left,border_left],[-50,50],'--k');
plot([0,0],[-50,50],'--k');
plot([border_right,border_right],[-50,50],'--k');
% drawrectangle('Position',[-5,-5,10,10]);

axis('equal');
for t = 1:T+1
    if linear_dyn
        ext_a = viscircles([X(1,t),X(2,t)],car_rad,'Color','b');
        ext_ap = viscircles([Xp(1,t),Xp(2,t)],car_rad,'Color','b','LineStyle',':');
        ext_b = viscircles([X(5,t),X(6,t)],car_rad,'Color','r');
        ext_c = viscircles([X(9,t),X(10,t)],car_rad,'Color','g');
        ext_cp = viscircles([Xp(9,t),Xp(10,t)],car_rad,'Color','g','LineStyle',':');
        ext_d = viscircles([X(13,t),X(14,t)],car_rad,'Color','c');
        ext_dp = viscircles([Xp(13,t),Xp(14,t)],car_rad,'Color','c','LineStyle',':');
        
%         col_a = viscircles([X(1,t),X(2,t)],min_dist/2,'Color','k','LineStyle',':');
%         col_ap = viscircles([Xp(1,t),Xp(2,t)],min_dist/2,'Color','k','LineStyle',':');
%         col_b = viscircles([X(5,t),X(6,t)],min_dist/2,'Color','k','LineStyle',':');
%         col_c = viscircles([X(9,t),X(10,t)],min_dist/2,'Color','k','LineStyle',':');
%         col_d = viscircles([X(13,t),X(14,t)],min_dist/2,'Color','k','LineStyle',':');
        
%         path_a = plot(X(1,t),X(2,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
%         path_b = plot(X(5,t),X(6,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'r');
%         path_c = plot(X(9,t),X(10,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'g');
%         path_d = plot(X(13,t),X(14,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'g');
    else
        path_a = quiver(X(1,t),X(2,t),X(3,t)*cos(X(4,t)),X(3,t)*sin(X(4,t)),0,'b','linewidth',3,'MaxHeadSize',1);
        path_b = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
        path_c = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
    end
    axis([X(1,t)-20 X(1,t)+20 X(2,t)-20 X(2,t)+20]);
    pause(0.2);
%     delete(path_a);
%     delete(path_b);
%     delete(path_c);
%     delete(path_d);
    delete(ext_a);
    delete(ext_ap);
    delete(ext_b);
    delete(ext_c);
    delete(ext_cp);
    delete(ext_d);
    delete(ext_dp);
    
%     delete(col_a);
%     delete(col_ap);
%     delete(col_b);
%     delete(col_c);
%     delete(col_d);
end



%%

function const = avoid_constraint(x)
    const = zeros(0,1);
end

function const = empty_constraint(~)
    const = zeros(0,1);
end

function val = empty_cost(~)
    val = 0;
end

function const = inequality_constraint_a(x)
    global border_left border_right;
    global sq_dist min_dist;
    xa = x(1);
    ya = x(2);
    va = x(4);
    
    const = [-xa + (border_left+min_dist/2);
             -(border_right-min_dist/2) + xa;
             -va];
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
    xd = x(13:14);
    const_ab = (xa-xb)'*(xa-xb)-sq_dist;
    const_ac = (xa-xc)'*(xa-xc)-sq_dist;
    const_ad = (xa-xd)'*(xa-xd)-sq_dist;
    const = [const;
             -const_ab];
end

function const = inequality_constraint_b(x)
    global border_left border_right;
    global sq_dist min_dist;
    xb = x(5);
    yb = x(6);
    const = [-xb + (border_left+min_dist/2);
             -(border_right-min_dist/2) + xb];
%     xa = x(1:2);
%     xb = x(5:6);
%     xc = x(9:10);
%     xd = x(13:14);
%     const_ab = (xa-xb)'*(xa-xb)-sq_dist;
%     const_bc = (xb-xc)'*(xb-xc)-sq_dist;
%     const_bd = (xb-xd)'*(xb-xd)-sq_dist;
%     const = [-const;
%              -const_bc;
%              -const_bd];
end

function const = inequality_constraint_c(x)
    global border_left border_right
    global sq_dist min_dist;
    xc = x(9);
    yc = x(10);
    const = [-xc + (border_left+min_dist/2);
             -(border_right-min_dist/2) + xc];
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
%     xd = x(13:14);
    const_ac = (xa-xc)'*(xa-xc)-sq_dist;
    const_bc = (xb-xc)'*(xb-xc)-sq_dist;
%     const_ad = (xa-xd)'*(xa-xd)-sq_dist;
    const = [const;
            -const_ac;
            -const_bc];
end

function const = inequality_constraint_d(x)
    global border_left border_right
    global sq_dist min_dist;
    xd = x(13);
    yd = x(14);
    const = [-xd + (border_left+min_dist/2);
             -(border_right-min_dist/2) + xd];
    xa = x(1:2);
    xb = x(5:6);
%     xc = x(9:10);
    xd = x(13:14);
    const_ad = (xa-xd)'*(xa-xd)-sq_dist;
    const_bd = (xb-xd)'*(xb-xd)-sq_dist;
    const = [const;
             -const_ad;
             -const_bd];
end

function const = terminal_constraint_a(x)
    const = zeros(0,1);
    global xga;
    xa = x(1:4);
    const = xa(1)-xga(1);
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
    const = xc(1)-xgc(1);
end

function const = terminal_constraint_d(x)
    const = zeros(0,1);
    global xgd;
    xd = x(13:16);
    const = xd(1)-xgd(1);
end

function const = terminal_ineq_constraint_a(x)
    const = zeros(0,1);
    xa = x(1:4);
    const = [xa(2)-25;
             xa(1)-2];
    const = [const; inequality_constraint_a(x)];
end

function const = mid_ineq_constraint_a(x)
    const = zeros(0,1);
    xa = x(1:4);
    const = [-xa(1)];
    const = [const; inequality_constraint_a(x)];
end

function const = terminal_ineq_constraint_b(x)
    const = zeros(0,1);
    xb = x(5:8);
    const = [-xb(1);
             xb(2)-20];
    const = [const; inequality_constraint_b(x)];
end

function const = terminal_ineq_constraint_c(x)
    const = zeros(0,1);
    xc = x(9:12);    
    const = [xc(2)-10];
    const = [const; inequality_constraint_c(x)];
end

function const = terminal_ineq_constraint_d(x)
    const = zeros(0,1);
    xd = x(13:16);
    const = [xd(2)-15];
    const = [const; inequality_constraint_d(x)];
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

function val = running_cost_a_polite(xu)
    global probc;
    ua = xu(17:18);
    uc = xu(21:22);
    yc = xu(10);
    yd = xu(14);
    val = 10*ua'*ua + 2*(probc*running_cost_c(xu) + (1-probc)*running_cost_d(xu));
end

function val = running_cost_a(xu)
    ua = xu(17:18);
    uc = xu(21:22);
    yc = xu(10);
    yd = xu(14);
    val = 10*ua'*ua;
end

function val = running_cost_b(xu)
    ub = xu(19:20);
    val = 10*ub'*ub;
end

function val = running_cost_c(xu)
    uc = xu(21:22);
    yc = xu(10);
    val = 1*uc'*uc-0.5*yc;
end

function val = running_cost_d(xu)
    ud = xu(23:24);
    yd = xu(14);
    val = 1*ud'*ud-1*yd;
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
    xd = xu(13:16);
    ua = xu(17:18);
    ub = xu(19:20);
    uc = xu(21:22);
    ud = xu(23:24);
    
    za = xa + 0.1*dubins_dyn(xa,ua);
    zb = xb + 0.1*dubins_dyn(xb,ub);
    zc = xc + 0.1*dubins_dyn(xc,uc);
    zd = xd + 0.1*dubins_dyn(xd,ud);
    % maybe wrap angles( if it doesn't break differentiability)
    next = [za;zb;zc;zd];
end