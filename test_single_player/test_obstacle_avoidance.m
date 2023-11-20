%% test_sqp_solver

clear; clc; close all;
global linear_dyn;
global shape width;
width = 0.75;
shape = 'square';
linear_dyn = true;
global xga xgb xgc;
global half_border;
half_border = 6;
global sq_dist;

min_dist = 1;
sq_dist = min_dist*min_dist;

m{1} = 2;
m{2} = 2;
m{3} = 2;
m{4} = 12;
N = 4;
n = 12;
T = 40;

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
    l{t,4} = @running_cost_r;
    
    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;
    h{t,4} = @empty_constraint;
    
    
    g{t,1} = @inequality_constraint_a;
    g{t,2} = @inequality_constraint_b;
    g{t,3} = @inequality_constraint_c;
    g{t,4} = @inequality_constraint_r;

end
l{T+1,1} = @terminal_cost_a;
l{T+1,2} = @terminal_cost_b;
l{T+1,3} = @terminal_cost_c;
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;
l{T+1,4} = @empty_cost;

h{T+1,1} = @terminal_constraint_a;
h{T+1,2} = @terminal_constraint_b;
h{T+1,3} = @terminal_constraint_c;
% h{T+1,1} = @empty_constraint;
h{T+1,4} = @empty_constraint;

g{T+1,1} = @inequality_constraint_a;
g{T+1,2} = @inequality_constraint_b;
g{T+1,1} = @empty_constraint;
g{T+1,2} = @empty_constraint;
g{T+1,3} = @empty_constraint;
g{T+1,4} = @empty_constraint;

rr{1} = false;
rr{2} = false;
rr{3} = false;
rr{4} = true;

[residuals, X,U] = sqp_solver(@f, @avoid_constraint, h, g, l, n, m, N, T, rr, x0);

disp('here');
%% plot solution
close all;
make_video = true;
if make_video
    writerObj = VideoWriter('avoid_ol','MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    axis tight
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    hold on;

else
    figure;
    hold on;
end

axis([-2*half_border 2*half_border -2*half_border 2*half_border]);
boundary = rectangle('Position',[-half_border,-half_border,2*half_border,2*half_border]);
axis('equal');
for t = 1:T+1
    if linear_dyn
%         ext_a = rectangle('Position',[X(1,t)-width,X(2,t)-width,2*width,2*width]);
%         ext_b = rectangle('Position',[X(5,t)-width,X(6,t)-width,2*width,2*width]);
%         ext_c = rectangle('Position',[X(9,t)-width,X(10,t)-width,2*width,2*width]);
        
        ext_a = viscircles([X(1,t),X(2,t)],min_dist/2,'Color','b');
        ext_b = viscircles([X(5,t),X(6,t)],min_dist/2,'Color','r');
        ext_c = viscircles([X(9,t),X(10,t)],min_dist/2,'Color','g');
        
%         if t < T+1
%             spot_ab = plot(U{4}(1,t),U{4}(2,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'g');
%             spot_ac = plot(U{4}(3,t),U{4}(4,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r');
%             spot_bc = plot(U{4}(5,t),U{4}(6,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'g');
%             spot_ba = plot(U{4}(7,t),U{4}(8,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b');
%             spot_cb = plot(U{4}(9,t),U{4}(10,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'r');
%             spot_ca = plot(U{4}(11,t),U{4}(12,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');
%         end
%         path_a = plot(X(1,t),X(2,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
%         path_b = plot(X(5,t),X(6,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'r');
%         path_c = plot(X(9,t),X(10,t),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'g');
    else
        path_a = quiver(X(1,t),X(2,t),X(3,t)*cos(X(4,t)),X(3,t)*sin(X(4,t)),0,'b','linewidth',3,'MaxHeadSize',1);
        path_b = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
        path_c = quiver(X(5,t),X(6,t),X(7,t)*cos(X(8,t)),X(7,t)*sin(X(8,t)),0,'r','linewidth',3,'MaxHeadSize',1);
    end
    axis([-10 10 -10 10]);
    if make_video
        frame = getframe;
        writeVideo(writerObj,frame);
    else
        pause(0.1);
    end
%     delete(path_a);
%     delete(path_b);
%     delete(path_c);
    delete(ext_a);
    delete(ext_b);
    delete(ext_c);
%     if t < T+1
%         delete(spot_ab);
%         delete(spot_ac);
%         delete(spot_bc);
%         delete(spot_ba);
%         delete(spot_cb);
%         delete(spot_ca);
%     end
end

if make_video
    close(writerObj);
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

function const = inequality_constraint_a(xu)
    global width;
    global half_border;
    
    xa = xu(1);
    ya = xu(2);
    const = [xa + half_border;
            -xa + half_border;
             ya + half_border;
            -ya + half_border];
    
    xa = xu(1:2);
    xb = xu(5:6);
    xc = xu(9:10);
    ur = xu(19:30);
    rad2 = 2*width^2;
% 
    dac = xa-ur(7:8);
    const = [const];
%     const_ab = (xa-xb)'*(xa-xb)-30;
    const_ac = (xa-xc)'*(xa-xc)-30;
%     const_bc = (xb-xc)'*(xb-xc)-30;
    const = [const;
            const_ac];
%              const_ac;
%              const_bc];
    
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

function const = inequality_constraint_r(xu)
    xa = xu(1:2);
    xb = xu(5:6);
    xc = xu(9:10);
    ur = xu(19:30);
    
    const = [point_in_box(ur(1:2),xa);
           point_in_box(ur(3:4),xb);
           point_in_box(ur(5:6),xa);
           point_in_box(ur(7:8),xc);
           point_in_box(ur(9:10),xb);
           point_in_box(ur(11:12),xc)];
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

function val = running_cost_r(xu)
    xa = xu(1:2);
    xb = xu(5:6);
    xc = xu(9:10);
    ur = xu(19:30);
    
    dab = xa-ur(3:4);
    dba = xb-ur(1:2);
    
    dac = xa-ur(7:8);
    dca = xc-ur(5:6);
    
    dbc = xb-ur(11:12);
    dcb = xc-ur(9:10);
    
    val = dab'*dab + dba'*dba + dac'*dac + dca'*dca + dbc'*dbc + dcb'*dcb;
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

function con = point_in_box(point, center)
    global width;
    con = [point(1) - (center(1)-width);
           center(1)+width-point(1);
           point(2) - (center(2)-width);
           center(2)+width-point(2)];
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