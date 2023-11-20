%% no end goal driving game

clear; clc; close all;

global dt;
dt = 0.1;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;

border_left = -3.7;
border_right = 3.7;


global avoid_r half_l half_w;
half_l = 1.5;
half_w = 0.75;
avoid_r = 2*sqrt(half_l^2+half_w^2);

global prob_c;
prob_c = .9;

global gv_a gv_b gv_c1 gv_c2;

gv_a = 5;
gv_b = 7.5;
gv_c1 = 5;
gv_c2 = 5;

x0 = [border_left/2;-5;0;gv_a;
       border_right/2;1;0;gv_a;
       border_right/2;10;0;gv_a;
       border_right/2;15;0;gv_a];

m{1} = 2;
m{2} = 2;
m{3} = 2;
m{4} = 2;
N = 4;
n = 16;
global T
T = 100;

for t = 1:T
    l{t,1} = @running_cost_a;
    l{t,2} = @running_cost_b;
    l{t,3} = @running_cost_c1;
    l{t,4} = @running_cost_c2;
  
    g{t,1} = @running_ineq_constraint_a;
    g{t,2} = @running_ineq_constraint_b;
    g{t,3} = @running_ineq_constraint_c1;
    g{t,4} = @running_ineq_constraint_c2;

    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;
    h{t,4} = @empty_constraint;

end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;
l{T+1,4} = @empty_cost;

g{T+1,1} = @running_ineq_constraint_a;
g{T+1,2} = @running_ineq_constraint_b;
g{T+1,3} = @running_ineq_constraint_c1;
g{T+1,4} = @running_ineq_constraint_c2;

% h{T+1,1} = @terminal_constraint_a;
% h{T+1,2} = @terminal_constraint_b;
% h{T+1,3} = @terminal_constraint_c1;
% h{T+1,4} = @terminal_constraint_c2;

h{T+1,1} = @empty_constraint;
h{T+1,2} = @empty_constraint;
h{T+1,3} = @empty_constraint;
h{T+1,4} = @empty_constraint;

params.tauval_init = {1,1,0.001,0.001};
params.tauval_decrease = 0.1;
params.tauval_tolerance = 0.001;
params.resid_tolerance = 0.001;
params.debug_plot = true;
params.max_linesearch_iters = 1;
params.alpha_init = 1;
params.alpha_init2 = 1;
params.beta = 0.5;
params.open_loop = true;
params.crash_start = false;
params.independent_updates = true;

[residuals, xval] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);

% params.open_loop = false;
% [residuals, xval_p] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);

%%
plot_solution(xval,x0);


%%

function plot_solution(xval,x0,varargin)
    make_vid = false;
    if length(varargin) > 0
        plot_polite = true;
    else
        plot_polite = false;
    end
    close all;
    fh = figure; 
    hold on;
    if make_vid
        writerObj = VideoWriter('possible_merge','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    regi = [-20,20,-20,20];
    global border_left border_right avoid_r linear_dyn T half_l half_w;
    
    plot([border_left, border_left],[-100,100],'-k');
    plot([border_left/2, border_left/2],[-100,100],'--k');
    plot([0, 0],[-100,100],'-k');
    plot([border_right/2, border_right/2],[-100,100],'--k');
    plot([border_right, border_right],[-100,100],'-k');

    xx(1,:) = x0;
    xxp(1,:) = x0;
    if linear_dyn
        spot_a = viscircles([xx(1,1),xx(1,2)],avoid_r/2,'Color', 'b');
        car_a = patch('XData',[xx(1,1)-half_w,xx(1,1)-half_w,xx(1,1)+half_w,xx(1,1)+half_w,xx(1,1)-half_w],...
                      'YData',[xx(1,2)-half_l,xx(1,2)+half_l,xx(1,2)+half_l,xx(1,2)-half_l,xx(1,2)-half_l],...
                      'FaceAlpha',1,'FaceColor', 'b');
        spot_b = viscircles([xx(1,5),xx(1,6)],avoid_r/2,'Color', 'r');
        car_b = patch('XData',[xx(1,5)-half_w,xx(1,5)-half_w,xx(1,5)+half_w,xx(1,5)+half_w,xx(1,5)-half_w],...
                       'YData',[xx(1,6)-half_l,xx(1,6)+half_l,xx(1,6)+half_l,xx(1,6)-half_l,xx(1,6)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'r');
        spot_c1 = viscircles([xx(1,9),xx(1,10)],avoid_r/2,'Color', 'g');
        car_c1 = patch('XData',[xx(1,9)-half_w,xx(1,9)-half_w,xx(1,9)+half_w,xx(1,9)+half_w,xx(1,9)-half_w],...
                       'YData',[xx(1,10)-half_l,xx(1,10)+half_l,xx(1,10)+half_l,xx(1,10)-half_l,xx(1,10)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'g');
        spot_c2 = viscircles([xx(1,13),xx(1,14)],avoid_r/2,'Color', 'c');
        car_c2 = patch('XData',[xx(1,13)-half_w,xx(1,13)-half_w,xx(1,13)+half_w,xx(1,13)+half_w,xx(1,13)-half_w],...
                       'YData',[xx(1,14)-half_l,xx(1,14)+half_l,xx(1,14)+half_l,xx(1,14)-half_l,xx(1,14)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'g');
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
    axis([xx(1,1)-20,xx(1,1)+20,xx(1,2)-20,xx(1,2)+20]);
    axis('equal');
    if make_vid
        frame = getframe(fh);
        writeVideo(writerObj,frame);
    else
        pause(0.01);
    end
    delete(spot_a);
    delete(spot_b);
    delete(spot_c1);
    delete(spot_c2);
    delete(car_a);
    delete(car_b);
    delete(car_c1);
    delete(car_c2);
    for t = 1:T
        xx(t+1,:) = xval{t+1};
        if plot_polite
            xxp(t+1,:) = varargin{1}{t+1};
        end
        if linear_dyn
            spot_a = viscircles([xx(t+1,1),xx(t+1,2)],avoid_r/2,'Color', 'b');
            car_a = patch('XData',[xx(t+1,1)-half_w,xx(t+1,1)-half_w,xx(t+1,1)+half_w,xx(t+1,1)+half_w,xx(t+1,1)-half_w],...
                   'YData',[xx(t+1,2)-half_l,xx(t+1,2)+half_l,xx(t+1,2)+half_l,xx(t+1,2)-half_l,xx(t+1,2)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'b');
            spot_b = viscircles([xx(t+1,5),xx(t+1,6)],avoid_r/2,'Color', 'r');
            car_b = patch('XData',[xx(t+1,5)-half_w,xx(t+1,5)-half_w,xx(t+1,5)+half_w,xx(t+1,5)+half_w,xx(t+1,5)-half_w],...
                   'YData',[xx(t+1,6)-half_l,xx(t+1,6)+half_l,xx(t+1,6)+half_l,xx(t+1,6)-half_l,xx(t+1,6)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'r');
            spot_c1 = viscircles([xx(t+1,9),xx(t+1,10)],avoid_r/2,'Color', 'g');
            car_c1 = patch('XData',[xx(t+1,9)-half_w,xx(t+1,9)-half_w,xx(t+1,9)+half_w,xx(t+1,9)+half_w,xx(t+1,9)-half_w],...
                       'YData',[xx(t+1,10)-half_l,xx(t+1,10)+half_l,xx(t+1,10)+half_l,xx(t+1,10)-half_l,xx(t+1,10)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'g');
            spot_c2 = viscircles([xx(t+1,13),xx(t+1,14)],avoid_r/2,'Color', 'c');
            car_c2 = patch('XData',[xx(t+1,13)-half_w,xx(t+1,13)-half_w,xx(t+1,13)+half_w,xx(t+1,13)+half_w,xx(t+1,13)-half_w],...
                       'YData',[xx(t+1,14)-half_l,xx(t+1,14)+half_l,xx(t+1,14)+half_l,xx(t+1,14)-half_l,xx(t+1,14)-half_l],...
                       'FaceAlpha',1,'FaceColor', 'g');
            if plot_polite
                pspot_a = viscircles([xxp(t+1,1),xxp(t+1,2)],avoid_r/2,'Color', 'b','LineStyle',':');
                pcar_a = patch('XData',[xxp(t+1,1)-half_w,xxp(t+1,1)-half_w,xxp(t+1,1)+half_w,xxp(t+1,1)+half_w,xxp(t+1,1)-half_w],...
                       'YData',[xxp(t+1,2)-half_l,xxp(t+1,2)+half_l,xxp(t+1,2)+half_l,xxp(t+1,2)-half_l,xxp(t+1,2)-half_l],...
                           'FaceAlpha',1,'FaceColor', 'b');
                pspot_b = viscircles([xxp(t+1,5),xxp(t+1,6)],avoid_r/2,'Color', 'r','LineStyle',':');
                pcar_b = patch('XData',[xxp(t+1,5)-half_w,xxp(t+1,5)-half_w,xxp(t+1,5)+half_w,xxp(t+1,5)+half_w,xxp(t+1,5)-half_w],...
                       'YData',[xxp(t+1,6)-half_l,xxp(t+1,6)+half_l,xxp(t+1,6)+half_l,xxp(t+1,6)-half_l,xxp(t+1,6)-half_l],...
                           'FaceAlpha',1,'FaceColor', 'r');
                pspot_c1 = viscircles([xxp(t+1,9),xxp(t+1,10)],avoid_r/2,'Color', 'g','LineStyle',':');
                pcar_c1 = patch('XData',[xxp(t+1,9)-half_w,xxp(t+1,9)-half_w,xxp(t+1,9)+half_w,xxp(t+1,9)+half_w,xxp(t+1,9)-half_w],...
                           'YData',[xxp(t+1,10)-half_l,xxp(t+1,10)+half_l,xxp(t+1,10)+half_l,xxp(t+1,10)-half_l,xxp(t+1,10)-half_l],...
                           'FaceAlpha',1,'FaceColor', 'g');
                pspot_c2 = viscircles([xxp(t+1,13),xxp(t+1,14)],avoid_r/2,'Color', 'c','LineStyle',':');
                pcar_c2 = patch('XData',[xxp(t+1,13)-half_w,xxp(t+1,13)-half_w,xxp(t+1,13)+half_w,xxp(t+1,13)+half_w,xxp(t+1,13)-half_w],...
                           'YData',[xxp(t+1,14)-half_l,xxp(t+1,14)+half_l,xxp(t+1,14)+half_l,xxp(t+1,14)-half_l,xxp(t+1,14)-half_l],...
                           'FaceAlpha',1,'FaceColor', 'g');

            end
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
        axis([xx(t+1,9)-20,xx(t+1,9)+20,xx(t+1,10)-20,xx(t+1,10)+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.01);
        end
        delete(spot_a);
        delete(spot_b);
        delete(spot_c1);
        delete(spot_c2);
        delete(car_a);
        delete(car_b);
        delete(car_c1);
        delete(car_c2);
        if plot_polite
            delete(pspot_a);
            delete(pspot_b);
            delete(pspot_c1);
            delete(pspot_c2);
            delete(pcar_a);
            delete(pcar_b);
            delete(pcar_c1);
            delete(pcar_c2);
        end
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

function val = empty_cost(~)
    val = 0;
end

function const = running_ineq_constraint_a(x)
    global border_left border_right half_w;
    global avoid_r;
    xa = x(1:2);
    xb = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    
    const = [(xa-xb)'*(xa-xb)-avoid_r*avoid_r;
             (xa-xc1)'*(xa-xc1)-avoid_r*avoid_r;
             (xa-xc2)'*(xa-xc2)-avoid_r*avoid_r];
    const = [const;
             border_right-half_w-x(1);
             x(1)-border_left-half_w];
end

function const = running_ineq_constraint_b(x)
    global border_left border_right half_w;
    global avoid_r;
    xa = x(1:2);
    xb = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    const = [(xb-xc1)'*(xb-xc1)-avoid_r*avoid_r;
             (xb-xc2)'*(xb-xc2)-avoid_r*avoid_r];
    const = [const;
             border_right-half_w-x(5);
             x(5)-border_left-half_w];
end

function const = running_ineq_constraint_c1(x)
    global border_left border_right half_w;
    global avoid_r;
    xa = x(1:2);
    xb = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    const = [(xc1-xc2)'*(xc1-xc2)-avoid_r*avoid_r];
    const = [const;
            border_right-half_w-x(9);
             x(9)-border_left-half_w];
end

function const = running_ineq_constraint_c2(x)
    global border_left border_right half_w;
    global avoid_r;
    xa = x(1:2);
    xb = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    const = zeros(0,1);
    const = [border_right-half_w-x(13);
             x(13)-border_left-half_w];
end

function const = terminal_constraint_a(x)
    global goal_a;
    xa = x(1:2);
    const = [xa-goal_a];
end

function const = terminal_constraint_b(x)
    global goal_b;
    xb = x(5:6);
    const = [xb-goal_b];
end

function const = terminal_constraint_c1(x)
    global goal_c1;
    xc1 = x(9:10);
    const = [xc1-goal_c1];
end

function const = terminal_constraint_c2(x)
    global goal_c2;
    xc2 = x(13:14);
    const = [xc2-goal_c2];
end

function val = terminal_cost(~)
    val = 0;
end

function val = running_cost_a(xu)
    global gv_a;
    xa = xu(1:4);
    ua = xu(end-8+1:end-8+2);
    ub = xu(end-6+1:end-6+2);
    uc1 = xu(end-4+1:end-4+2);
    uc2 = xu(end-2+1:end-2+2);
    val = 1*ua'*ua + 1*(xa(4)-gv_a)^2+lane_cost(xa);
end

function val = running_cost_a_polite(xu)
    global prob_c;
    val = running_cost_a(xu) + prob_c*running_cost_c1(xu) + (1-prob_c)*running_cost_c2(xu);
end

function val = running_cost_b(xu)
    global gv_b;
    xb = xu(5:8);
    ua = xu(end-8+1:end-8+2);
    ub = xu(end-6+1:end-6+2);
    uc1 = xu(end-4+1:end-4+2);
    uc2 = xu(end-2+1:end-2+2);
    val = 1*ub'*ub + 10*(xb(4)-gv_b)^2+lane_cost(xb);
end

function val = running_cost_c1(xu)
    global gv_c1;
    xc1 = xu(9:12);
    ua = xu(end-8+1:end-8+2);
    ub = xu(end-6+1:end-6+2);
    uc1 = xu(end-4+1:end-4+2);
    uc2 = xu(end-2+1:end-2+2);
    val = 1*uc1'*uc1 + 1*(xc1(4)-gv_c1)^2+lane_cost(xc1);
end

function val = running_cost_c2(xu)
    global gv_c2;
    xc2 = xu(13:16);
    ua = xu(end-8+1:end-8+2);
    ub = xu(end-6+1:end-6+2);
    uc1 = xu(end-4+1:end-4+2);
    uc2 = xu(end-2+1:end-2+2);
    val = 1*uc2'*uc2 + 1*(xc2(4)-gv_c2)^2+lane_cost(xc2);
end

function val = running_cost_c1_polite(xu)
    val = running_cost_c1(xu) + running_cost_a(xu);
end

function val = running_cost_c2_polite(xu)
    val = running_cost_c1(xu) + running_cost_a(xu);
end

function val = lane_cost(x)
    global border_left border_right;
    pos = x(1);
    val_left = (border_left-pos)^2;
    val_right= (pos-border_right)^2;
    val1 = if_else(pos<border_left,val_left,0);
    val2 = if_else(border_right<pos,val_right,0);
    val = 0;
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
    xc1 = xu(9:12);
    xc2 = xu(13:16);
    ua = xu(17:18);
    ub = xu(19:20);
    uc1 = xu(21:22);
    uc2 = xu(23:24);
    
    za = xa + dt*unicycle_dyn(xa,ua);
    zb = xb + dt*unicycle_dyn(xb,ub);
    zc1 = xc1 + dt*unicycle_dyn(xc1,uc1);
    zc2 = xc2 + dt*unicycle_dyn(xc2,uc2);
    next = [za;zb;zc1;zc2];
end