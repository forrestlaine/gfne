%% game of clones

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

global prob_a1 prob_b1 prob_c1;
prob_a1 = 0.9;
prob_b1 = 0.1;
prob_c1 = 0.9;

global gv_a1 gv_a2 gv_b1 gv_b2 gv_c1 gv_c2;

gv_a1 = 5;
gv_a2 = 5;
gv_b1 = 3;
gv_b2 = 5;
gv_c1 = 7;
gv_c2 = 7;

x0 = [border_right/2;0;0;gv_a1;
      border_right/2;0;0;gv_a2;
       border_right/2;10;0;gv_a1;
       border_right/2;10;0;gv_a2;
       border_left/2;-10;0;gv_a1;
       border_left/2;-10;0;gv_a2];

m{1} = 2;
m{2} = 2;
m{3} = 2;
m{4} = 2;
m{5} = 2;
m{6} = 2;
N = 6;
n = N*4;
global T
T = 50;

for t = 1:T
    l{t,1} = @running_cost_a1_polite;
    l{t,2} = @running_cost_a2_polite;
    l{t,3} = @running_cost_b1_polite;
    l{t,4} = @running_cost_b2_polite;
    l{t,5} = @running_cost_c1_polite;
    l{t,6} = @running_cost_c2_polite;
  
    g{t,1} = @running_ineq_constraint_a1;
    g{t,2} = @running_ineq_constraint_a2;
    g{t,3} = @running_ineq_constraint_b1;
    g{t,4} = @running_ineq_constraint_b2;
    g{t,5} = @running_ineq_constraint_c1;
    g{t,6} = @running_ineq_constraint_c2;

    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;
    h{t,4} = @empty_constraint;
    h{t,5} = @empty_constraint;
    h{t,6} = @empty_constraint;

end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;
l{T+1,4} = @empty_cost;
l{T+1,5} = @empty_cost;
l{T+1,6} = @empty_cost;

g{T+1,1} = @terminal_ineq_constraint_a1;
g{T+1,2} = @terminal_ineq_constraint_a2;
g{T+1,3} = @terminal_ineq_constraint_b1;
g{T+1,4} = @terminal_ineq_constraint_b2;
g{T+1,5} = @terminal_ineq_constraint_c1;
g{T+1,6} = @terminal_ineq_constraint_c2;

h{T+1,1} = @terminal_constraint_a1;
h{T+1,2} = @terminal_constraint_a2;
h{T+1,3} = @terminal_constraint_b1;
h{T+1,4} = @terminal_constraint_b2;
h{T+1,5} = @terminal_constraint_c1;
h{T+1,6} = @terminal_constraint_c2;


params.tauval_init = {.1,.1,.1,.1,.1,.1};
params.tauval_decrease = 0.1;
params.tauval_tolerance = 0.001;
params.resid_tolerance = 0.001;
params.debug_plot = true;
params.max_linesearch_iters = 2;
params.max_cl_linesearch_iters = 3;
params.alpha_init = 1;
params.beta = 0.5;
params.open_loop = true;
params.crash_start = false;
params.independent_updates = true;

[residuals, xval] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);

% params.open_loop = true;
% [residuals, xval_ol] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);

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

    spot_a1 = viscircles(xx(1,1:2),avoid_r/2,'Color', 'b');
    car_a1 = draw_car(xx(1,1:2),'b');
    spot_a2 = viscircles(xx(1,5:6),avoid_r/2,'Color', 'b');
    car_a2 = draw_car(xx(1,5:6),'b');

    spot_b1 = viscircles(xx(1,9:10),avoid_r/2,'Color', 'r');
    car_b1 = draw_car(xx(1,9:10),'r');
    spot_b2 = viscircles(xx(1,13:14),avoid_r/2,'Color', 'r');
    car_b2 = draw_car(xx(1,13:14),'r');

    spot_c1 = viscircles(xx(1,17:18),avoid_r/2,'Color', 'g');
    car_c1 = draw_car(xx(1,17:18),'g');
    spot_c2 = viscircles(xx(1,21:22),avoid_r/2,'Color', 'g');
    car_c2 = draw_car(xx(1,21:22),'g');
    
    axis([xx(1,9)-20,xx(1,9)+20,xx(1,10)-20,xx(1,10)+20]);
    axis('equal');
    if make_vid
        frame = getframe(fh);
        writeVideo(writerObj,frame);
    else
        pause(0.01);
    end
    delete(spot_a1);
    delete(car_a1);
    delete(spot_a2);
    delete(car_a2);
    delete(spot_b1);
    delete(car_b1);
    delete(spot_b2);
    delete(car_b2);
    delete(spot_c1);
    delete(car_c1);
    delete(spot_c2);
    delete(car_c2);
    
    for t = 1:T
        xx(t+1,:) = xval{t+1};
        if plot_polite
            xxp(t+1,:) = varargin{1}{t+1};
        end
        
        spot_a1 = viscircles(xx(t+1,1:2),avoid_r/2,'Color', 'b');
        car_a1 = draw_car(xx(t+1,1:2),'b');
        spot_a2 = viscircles(xx(t+1,5:6),avoid_r/2,'Color', 'b');
        car_a2 = draw_car(xx(t+1,5:6),'b');

        spot_b1 = viscircles(xx(t+1,9:10),avoid_r/2,'Color', 'r');
        car_b1 = draw_car(xx(t+1,9:10),'r');
        spot_b2 = viscircles(xx(t+1,13:14),avoid_r/2,'Color', 'r');
        car_b2 = draw_car(xx(t+1,13:14),'r');

        spot_c1 = viscircles(xx(t+1,17:18),avoid_r/2,'Color', 'g');
        car_c1 = draw_car(xx(t+1,17:18),'g');
        spot_c2 = viscircles(xx(t+1,21:22),avoid_r/2,'Color', 'g');
        car_c2 = draw_car(xx(t+1,21:22),'g');

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
        
        axis([xx(t+1,9)-20,xx(t+1,9)+20,xx(t+1,10)-20,xx(t+1,10)+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.01);
        end
        delete(spot_a1);
        delete(car_a1);
        delete(spot_a2);
        delete(car_a2);
        delete(spot_b1);
        delete(car_b1);
        delete(spot_b2);
        delete(car_b2);
        delete(spot_c1);
        delete(car_c1);
        delete(spot_c2);
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

%% Running Constraints

function const = running_ineq_constraint_a1(x)
    global border_left border_right half_w;
    global avoid_r;
    
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [(xa1-xb1)'*(xa1-xb1)-avoid_r*avoid_r;
             (xa1-xb2)'*(xa1-xb2)-avoid_r*avoid_r];
    const = [const;
             border_right-half_w-xa1(1);
             xa1(1)-border_left-half_w];
end

function const = running_ineq_constraint_a2(x)
    global border_left border_right half_w;
    global avoid_r;
    
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [(xa2-xb1)'*(xa2-xb1)-avoid_r*avoid_r;
             (xa2-xb2)'*(xa2-xb2)-avoid_r*avoid_r];
    const = [const;
             border_right-half_w-xa2(1);
             xa2(1)-border_left-half_w];
end

function const = running_ineq_constraint_b1(x)
    global border_left border_right half_w;
    global avoid_r;
    
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [border_right-half_w-xb1(1);
             xb1(1)-border_left-half_w];
end

function const = running_ineq_constraint_b2(x)
    global border_left border_right half_w;
    global avoid_r;
    
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);

    const = [border_right-half_w-xb2(1);
             xb2(1)-border_left-half_w];
end


function const = running_ineq_constraint_c1(x)
    global border_left border_right half_w;
    global avoid_r;
    
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [(xc1-xa1)'*(xc1-xa1)-avoid_r*avoid_r;
             (xc1-xa2)'*(xc1-xa2)-avoid_r*avoid_r;
             (xc1-xb1)'*(xc1-xb1)-avoid_r*avoid_r;
             (xc1-xb2)'*(xc1-xb2)-avoid_r*avoid_r];
    const = [const;
             border_right-half_w-xc1(1);
             xc1(1)-border_left-half_w];
end

function const = running_ineq_constraint_c2(x)
    global border_left border_right half_w;
    global avoid_r;
    
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [(xc2-xa1)'*(xc2-xa1)-avoid_r*avoid_r;
             (xc2-xa2)'*(xc2-xa2)-avoid_r*avoid_r;
             (xc2-xb1)'*(xc2-xb1)-avoid_r*avoid_r;
             (xc2-xb2)'*(xc2-xb2)-avoid_r*avoid_r];
    const = [const;
             border_right-half_w-xc2(1);
             xc2(1)-border_left-half_w];
end

%% Terminal ineq Constraints

function const = terminal_ineq_constraint_a1(x)
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [running_ineq_constraint_a1(x);
             xa1(2)-xb1(2)];
end

function const = terminal_ineq_constraint_a2(x)
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [running_ineq_constraint_a2(x);
             xa2(2)-xb2(2)];
end

function const = terminal_ineq_constraint_b1(x)
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [running_ineq_constraint_b1(x)];
end

function const = terminal_ineq_constraint_b2(x)
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [running_ineq_constraint_b2(x)];
end


function const = terminal_ineq_constraint_c1(x)
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [running_ineq_constraint_c1(x)];
end

function const = terminal_ineq_constraint_c2(x)
    xa1 = x(1:2);
    xa2 = x(5:6);
    xb1 = x(9:10);
    xb2 = x(13:14);
    xc1 = x(17:18);
    xc2 = x(21:22);
    
    const = [running_ineq_constraint_c2(x)];
end

%% Terminal constraints 

function const = terminal_constraint_a1(x)
    xa1 = x(1:2);
    const = left_lane_constraint(xa1);
end

function const = terminal_constraint_a2(x)
    xa2 = x(5:6);
    const = left_lane_constraint(xa2);
end

function const = terminal_constraint_b1(x)
    xb1 = x(9:10);
    const = right_lane_constraint(xb1);
end

function const = terminal_constraint_b2(x)
    xb2 = x(13:14);
    const = right_lane_constraint(xb2);
end

function const = terminal_constraint_c1(x)
    xc1 = x(17:18);
    const = left_lane_constraint(xc1);
end

function const = terminal_constraint_c2(x)
    xc2 = x(21:22);
    const = right_lane_constraint(xc2);
end

function val = terminal_cost(~)
    val = 0;
end

%% Running costs
function val = running_cost_a1(xu)
    global gv_a1;
    xa1 = xu(1:4);
    ua1 = xu(25:26);
    val = ua1'*ua1 + 3*speed_cost(xa1,gv_a1);
end

function val = running_cost_a2(xu)
    global gv_a2;
    xa2 = xu(5:8);
    ua2 = xu(27:28);
    val = ua2'*ua2 + 3*speed_cost(xa2,gv_a2);
end

function val = running_cost_b1(xu)
    global gv_b1;
    xb1 = xu(9:12);
    ub1 = xu(29:30);
    val = ub1'*ub1 + 3*speed_cost(xb1,gv_b1);
end

function val = running_cost_b2(xu)
    global gv_b2;
    xb2 = xu(13:16);
    ub2 = xu(31:32);
    val = ub2'*ub2 + 3*speed_cost(xb2,gv_b2);
end

function val = running_cost_c1(xu)
    global gv_c1;
    xc1 = xu(17:20);
    uc1 = xu(33:34);
    val = uc1'*uc1 + 3*speed_cost(xc1,gv_c1);
end

function val = running_cost_c2(xu)
    global gv_c2;
    xc2 = xu(21:24);
    uc2 = xu(35:36);
    val = uc2'*uc2 + 3*speed_cost(xc2,gv_c2);
end

%% Polite costs
function val = running_cost_a1_polite(xu)
    global prob_c1;
    val = running_cost_a1(xu) + ...
        prob_c1*running_cost_c1(xu) + ...
        (1-prob_c1)*running_cost_c2(xu);
end

function val = running_cost_a2_polite(xu)
    global prob_c1;
    val = running_cost_a2(xu) + ...
        prob_c1*running_cost_c1(xu) + ...
        (1-prob_c1)*running_cost_c2(xu);
end

function val = running_cost_b1_polite(xu)
    global prob_a1;
    val = running_cost_b1(xu) + ...
        prob_a1*running_cost_a1(xu) + ...
        (1-prob_a1)*running_cost_a2(xu);
end

function val = running_cost_b2_polite(xu)
    global prob_a1;
    val = running_cost_b2(xu) + ...
        prob_a1*running_cost_a1(xu) + ...
        (1-prob_a1)*running_cost_a2(xu);
end

function val = running_cost_c1_polite(xu)
    val = running_cost_c1(xu);
end

function val = running_cost_c2_polite(xu)
    val = running_cost_c2(xu);
end
%% Helpers

function val = speed_cost(x,goal)
    val = (x(4)-goal)^2;
end

function const = left_lane_constraint(x)
    global border_left;
    const = x(1)-border_left/2;
end

function const = right_lane_constraint(x)
    global border_right;
    const = x(1)-border_right/2;
end

function pt = draw_car(center,color)
    global half_l half_w;
    pt = patch('XData',[center(1)-half_w,center(1)-half_w,center(1)+half_w,center(1)+half_w,center(1)-half_w],...
               'YData',[center(2)-half_l,center(2)+half_l,center(2)+half_l,center(2)-half_l,center(2)-half_l],...
               'FaceAlpha',1,'FaceColor', color);
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
    xa1 = xu(1:4);
    xa2 = xu(5:8);
    xb1 = xu(9:12);
    xb2 = xu(13:16);
    xc1 = xu(17:20);
    xc2 = xu(21:24);
    ua1 = xu(25:26);
    ua2 = xu(27:28);
    ub1 = xu(29:30);
    ub2 = xu(31:32);
    uc1 = xu(33:34);
    uc2 = xu(35:36);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    za2 = xa2 + dt*unicycle_dyn(xa2,ua2);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    zb2 = xb2 + dt*unicycle_dyn(xb2,ub2);
    zc1 = xc1 + dt*unicycle_dyn(xc1,uc1);
    zc2 = xc2 + dt*unicycle_dyn(xc2,uc2);
    next = [za1;za2;zb1;zb2;zc1;zc2];
end