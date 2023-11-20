%% test_ec_solver

clear; clc; close all;

global dt;
dt = 0.1;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;
global border_up;
global border_down;

border_left = -5;
border_right = 5;
border_up = 10;
border_down = -10;

global avoid_r;
avoid_r = 2;

global obstacle;
obstacle = [0;0];

global goal_a goal_b goal_c;
goal_a = [3;9];
goal_b = [-3;9];
goal_c = [-1;5];

m{1} = 2;
m{2} = 2;
m{3} = 2;
N = 3;
n = 12;
global T
T = 100;

x0 = [-2.5;-8;0;0;
       1;-7;0;0;
       -3;4;0;0];

l = cell(T+1,2);
h = cell(T+1,2);

for t = 1:T
    l{t,1} = @running_cost_a;
    l{t,2} = @running_cost_b;
    l{t,3} = @running_cost_c;
  
    g{t,1} = @empty_constraint;
    g{t,2} = @empty_constraint;
    g{t,3} = @empty_constraint; 

    g{t,1} = @running_ineq_constraint_a;
    g{t,2} = @running_ineq_constraint_b;
    g{t,3} = @running_ineq_constraint_c;

    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;

end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;

% l{T+1,1} = @terminal_cost_a;
% l{T+1,2} = @terminal_cost_b;
% l{T+1,3} = @terminal_cost_c;

% h{T+1,1} = @terminal_constraint_a;
% h{T+1,2} = @terminal_constraint_b;
% h{T+1,3} = @terminal_constraint_c;

% g{T+1,1} = @empty_constraint;
% g{T+1,2} = @empty_constraint;
% g{T+1,3} = @empty_constraint;

g{T+1,1} = @running_ineq_constraint_a;
g{T+1,2} = @running_ineq_constraint_b;
g{T+1,3} = @running_ineq_constraint_c;

h{T+1,1} = @terminal_ineq_constraint_a;
h{T+1,2} = @terminal_ineq_constraint_b;
h{T+1,3} = @terminal_ineq_constraint_c;

params.tauval_init = 1e-1;
params.tauval_decrease = 1e-2;
params.tauval_tolerance = 1e-3;
params.resid_tolerance = 1e-3;
params.debug_plot = false;
params.max_linesearch_iters = 1;
params.alpha_init = 1;
params.beta = 0.5;
params.open_loop = false;
params.crash_start = true;

[residuals, xval] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);
for t = 1:T
    l{t,2} = @running_cost_b_polite;
    l{t,3} = @running_cost_c_polite;
end
[residuals, xval_p] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);

%%
function plot_solution(xval,x0)
    global T goal_a goal_b goal_c half_l half_w
    make_vid = false;
    close all;
    fh = figure; 
    hold on;
    if make_vid
        writerObj = VideoWriter('ip_fb','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    viscircles(goal_a',1,'Color', 'b','LineStyle',':');
    viscircles(goal_b',1,'Color', 'r','LineStyle',':');
    viscircles(goal_c',1,'Color', 'g','LineStyle',':');
    viscircles([7.5, 5],1,'Color', 'k','LineStyle',':');
    viscircles([7.5, 0],1,'Color', 'k','LineStyle','-');
    viscircles([7.5, -5],1,'Color', 'k','LineStyle','-.');

    text(7.5,5,'Goal','HorizontalAlignment', 'center');
    text(7.5,0,'Selfish','HorizontalAlignment', 'center');
    text(7.5,-5,'Polite','HorizontalAlignment', 'center');

    global border_left border_down border_right border_up avoid_r linear_dyn;
    % viscircles(obstacle',1,'Color', 'k','LineStyle','-');
    rectangle('Position',[-10,-10,20,20]);
    rectangle('Position',[border_left, border_down, border_right-border_left, border_up-border_down]);
    half_w = .25;
    half_l = 0.5;

    xx(1,:) = x0;
    xxp(1,:) = x0;
    if linear_dyn
        spot_a = viscircles([xx(1,1),xx(1,2)],avoid_r/2,'Color', 'b');
        spot_b = viscircles([xx(1,5),xx(1,6)],avoid_r/2,'Color', 'r');
        spot_c = viscircles([xx(1,9),xx(1,10)],avoid_r/2,'Color', 'g');
    else
        spot_a = patch('XData',[xx(1,1)-half_w,xx(1,1)-half_w,xx(1,1)+half_w,xx(1,1)+half_w,xx(1,1)-half_w],...
                       'YData',[xx(1,2)-half_l,xx(1,2)+half_l,xx(1,2)+half_l,xx(1,2)-half_l,xx(1,2)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'b');
        spot_b = patch('XData',[xx(1,5)-half_w,xx(1,5)-half_w,xx(1,5)+half_w,xx(1,5)+half_w,xx(1,5)-half_w],...
                       'YData',[xx(1,6)-half_l,xx(1,6)+half_l,xx(1,6)+half_l,xx(1,6)-half_l,xx(1,6)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'r');
        spot_c = patch('XData',[xx(1,9)-half_w,xx(1,9)-half_w,xx(1,9)+half_w,xx(1,9)+half_w,xx(1,9)-half_w],...
                       'YData',[xx(1,10)-half_l,xx(1,10)+half_l,xx(1,10)+half_l,xx(1,10)-half_l,xx(1,10)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'g');

        rotate(spot_a, [0 0 1], rad2deg(xx(1,4))+90, [xx(1,1) xx(1,2) 0]);
        rotate(spot_b, [0 0 1], rad2deg(xx(1,8))+90, [xx(1,5) xx(1,6) 0]);
        rotate(spot_c, [0 0 1], rad2deg(xx(1,12))+90, [xx(1,9) xx(1,10) 0]);
    end

    %     spot_a = plot(xx(1,1), xx(1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    %     spot_b = plot(xx(1,5), xx(1,6),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    axis([-10 10 -10 10]);
    axis('equal');
    if make_vid
        frame = getframe(fh);
        writeVideo(writerObj,frame);
    else
        pause(0.025);
    end
    delete(spot_a);
    delete(spot_b);
    delete(spot_c);
    for t = 1:T
        xx(t+1,:) = xval{t+1};
%         xxp(t+1,:) = xval_p{t+1};
        if linear_dyn
            spot_a = viscircles([xx(t+1,1),xx(t+1,2)],avoid_r/2,'Color', 'b');
            spot_b = viscircles([xx(t+1,5),xx(t+1,6)],avoid_r/2,'Color', 'r');
            spot_c = viscircles([xx(t+1,9),xx(t+1,10)],avoid_r/2,'Color', 'g');
%             pspot_a = viscircles([xxp(t+1,1),xxp(t+1,2)],avoid_r/2,'Color', 'b','LineStyle','-.');
%             pspot_b = viscircles([xxp(t+1,5),xxp(t+1,6)],avoid_r/2,'Color', 'r','LineStyle','-.');
%             pspot_c = viscircles([xxp(t+1,9),xxp(t+1,10)],avoid_r/2,'Color', 'g','LineStyle','-.');
        else
            spot_a = patch('XData',[xx(t+1,1)-half_w,xx(t+1,1)-half_w,xx(t+1,1)+half_w,xx(t+1,1)+half_w,xx(t+1,1)-half_w],...
                   'YData',[xx(t+1,2)-half_l,xx(t+1,2)+half_l,xx(t+1,2)+half_l,xx(t+1,2)-half_l,xx(t+1,2)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'b');
            spot_b = patch('XData',[xx(t+1,5)-half_w,xx(t+1,5)-half_w,xx(t+1,5)+half_w,xx(t+1,5)+half_w,xx(t+1,5)-half_w],...
                   'YData',[xx(t+1,6)-half_l,xx(t+1,6)+half_l,xx(t+1,6)+half_l,xx(t+1,6)-half_l,xx(t+1,6)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'r');
            spot_c = patch('XData',[xx(t+1,9)-half_w,xx(t+1,9)-half_w,xx(t+1,9)+half_w,xx(t+1,9)+half_w,xx(t+1,9)-half_w],...
                       'YData',[xx(t+1,10)-half_l,xx(t+1,10)+half_l,xx(t+1,10)+half_l,xx(t+1,10)-half_l,xx(t+1,10)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'g');
            rotate(spot_a, [0 0 1], rad2deg(xx(t+1,4))+90, [xx(t+1,1) xx(t+1,2) 0]);
            rotate(spot_b, [0 0 1], rad2deg(xx(t+1,8))+90, [xx(t+1,5) xx(t+1,6) 0]);
            rotate(spot_c, [0 0 1], rad2deg(xx(t+1,12))+90, [xx(t+1,9) xx(t+1,10) 0]);
        end
    %         spot_a = plot(xx(t+1,1), xx(t+1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    %         spot_b = plot(xx(t+1,5), xx(t+1,6),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        axis([-10 10 -10 10]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.025);
        end
        delete(spot_a);
        delete(spot_b);
        delete(spot_c);
%         delete(pspot_a);
%         delete(pspot_b);
%         delete(pspot_c);
    end
    disp('plotted');
    if make_vid
        close(writerObj);
    end
end
%%

function const = empty_constraint(~)
    const = zeros(0,1);
end

function const = running_constraint_a(xu)
    va = xu(3);
    const = va-4;
end

function const = running_constraint_b(xu)
    vb = xu(7);
    const = vb-4;
end

function val = empty_cost(~)
    val = 0;
end

function const = running_ineq_constraint_a(x)
    global border_left border_right border_up border_down;
    global avoid_r obstacle goal_b goal_c;
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
%     const = [(xa-obstacle)'*(xa-obstacle)-avoid_r*avoid_r;
%             (xa-goal_b)'*(xa-goal_b)-avoid_r*avoid_r;
%             (xa-goal_c)'*(xa-goal_c)-avoid_r*avoid_r];
    const = [(xa-xb)'*(xa-xb)-avoid_r*avoid_r;
             (xa-xc)'*(xa-xc)-avoid_r*avoid_r];
    const = [const;
            border_right-x(1);
             x(1)-border_left;
             border_up-x(2);
             x(2)-border_down];
end

function const = running_ineq_constraint_b(x)
    global border_left border_right border_up border_down;
    global obstacle avoid_r goal_a goal_c;
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
%     const = [(xb-obstacle)'*(xb-obstacle)-avoid_r*avoid_r;
%              (xb-goal_a)'*(xb-goal_a)-avoid_r*avoid_r;
%              (xb-goal_c)'*(xb-goal_c)-avoid_r*avoid_r];
    const = (xc-xb)'*(xc-xb)-avoid_r*avoid_r;
    const = [const;
            border_right-x(5);
             x(5)-border_left;
             border_up-x(6);
             x(6)-border_down];
end

function const = running_ineq_constraint_c(x)
    global border_left border_right border_up border_down;
    global obstacle avoid_r goal_a goal_b;
    xa = x(1:2);
    xb = x(5:6);
    xc = x(9:10);
%     const = [(xc-obstacle)'*(xc-obstacle)-avoid_r*avoid_r;
%              (xc-goal_a)'*(xc-goal_a)-avoid_r*avoid_r;
%              (xc-goal_b)'*(xc-goal_b)-avoid_r*avoid_r];
    const = (xc-xa)'*(xc-xa)-avoid_r*avoid_r;
    const = [border_right-x(9);
             x(9)-border_left;
             border_up-x(10);
             x(10)-border_down];
end

function const = terminal_ineq_constraint_a(x)
    global goal_a;
    xa = x(1:2);
    const = [xa-goal_a; x(3:4)];
%     const = 1-(xa-goal)'*(xa-goal);
%     const = [x(2)-5;
%              x(1)];
end

function const = terminal_ineq_constraint_b(x)
    global goal_b;
    xb = x(5:6);
    const = [xb-goal_b; x(7:8)];
%     const = 1-(xb-goal)'*(xb-goal);
%     const = [x(6)-5;
%              x(5)-1];
end

function const = terminal_ineq_constraint_c(x)
    global goal_c;
    xc = x(9:10);
    const = [xc-goal_c; x(11:12)];
%     const = 1-(xc-goal)'*(xc-goal);
%     const = [x(10)-5;
%              -1-x(9)];
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
    const = zeros(0,1);
    xa = x(1:4);
    xb = x(5:8);
    const = xa-xb;
end

function const = terminal_constraint_b(x)
    const = zeros(0,1);
    xa = x(1:4);
    xb = x(5:8);
    xc = x(9:12);
    const = xb-xc;
%     const = const(1:4);
end

function const = terminal_constraint_c(x)
    const = zeros(0,1);
    goal = [-3;7;3;pi/2];
    xa = x(1:4);
    xb = x(5:8);
    xc = x(9:12);
    const = xa-goal;
end

function val = terminal_cost(~)
    val = 0;
end

function val = terminal_cost_a(x)
    goal = [-3;7;3;pi/2];
    xa = x(1:4);
    xb = x(5:8);
    xc = x(9:12);
    val = 100*(xb-goal)'*(xb-goal);
end

function val = terminal_cost_b(x)
    goal = [-3;7;3;pi/2];
    xa = x(1:4);
    xb = x(5:8);
    xc = x(9:12);
    val = 100*(xc-goal)'*(xc-goal);
end

function val = terminal_cost_c(x)
    goal = [-3;7;3;pi/2];
    xa = x(1:4);
    xb = x(5:8);
    xc = x(9:12);
    val = 100*(xa-goal)'*(xa-goal);
end

function val = running_cost_a(xu)
    ua = xu(end-6+1:end-6+2);
    val = 1*ua'*ua;
end

function val = running_cost_b(xu)
    ua = xu(end-6+1:end-6+2);
    ub = xu(end-4+1:end-4+2);
    val = 1*ub'*ub;
end

function val = running_cost_b_polite(xu)
    ua = xu(end-6+1:end-6+2);
    ub = xu(end-4+1:end-4+2);
    val = 1*ub'*ub + 1*ua'*ua;
end

function val = running_cost_c_polite(xu)
    ua = xu(end-6+1:end-6+2);
    ub = xu(end-4+1:end-4+2);
    uc = xu(end-2+1:end);
    val = 1*uc'*uc + 1*ub'*ub;
end

function val = running_cost_c(xu)
    uc = xu(end-2+1:end);
    val = 1*uc'*uc;
end

function vec = unicycle_dyn(x,u)
    global linear_dyn;
    vec = x(3)*cos(x(4));
    vec = [vec; x(3)*sin(x(4))];
    vec = [vec; u(1)];
    vec = [vec; u(2)];
    if linear_dyn
        vec = x(3);
        vec = [vec; x(4)];
        vec = [vec; u(1)];
        vec = [vec; u(2)];
    end
end

function next = f(xu)
    global dt;
    xa = xu(1:4);
    xb = xu(5:8);
    xc = xu(9:12);
    ua = xu(13:14);
    ub = xu(15:16);
    uc = xu(17:18);
    
    za = xa + dt*unicycle_dyn(xa,ua);
    zb = xb + dt*unicycle_dyn(xb,ub);
    zc = xc + dt*unicycle_dyn(xc,uc);
    % maybe wrap angles( if it doesn't break differentiability)
    %     za(4) = mod(za(4)+pi,2*pi)-pi;
    %     zb(4) = mod(zb(4)+pi,2*pi)-pi;
    next = [za;zb;zc];
end