%% merge in front of unknown speed
clear; clc; close all;

global dt;
dt = 0.25;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;

border_left = -4;
border_right = 4;

global center_lane_pen speed_pen;
center_lane_pen = 0;
speed_pen = 3;

global avoid_r half_l half_w;
half_l = 1.5;
half_w = 0.75;
avoid_r = 2*sqrt(half_l^2+half_w^2);

global prob_c1;

prob_c1 = 0.33;

global gv_a1 gv_b1 gv_c1 gv_c2;

gv_a1 = 5;
gv_b1 = 5;
gv_c1 = 5;
gv_c2 = 6;

x0 = [border_right/2; 0;0;gv_a1;
       border_right/2;5;0;gv_b1;
       border_left/2;-10;0;gv_c1;
       border_left/2;-10;0;gv_c2];

m{1} = 2;
m{2} = 2;
m{3} = 2;
m{4} = 2;
N = 4;
n = N*4;
global T
T = 40;

for t = 1:T
    l{t,1} = @running_cost_a1_polite;
    l{t,2} = @running_cost_b1_polite;
    l{t,3} = @running_cost_c1_polite;
    l{t,4} = @running_cost_c2_polite;
  
    g{t,1} = @running_ineq_constraint_a1;
    g{t,2} = @running_ineq_constraint_b1;
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

g{T+1,1} = @terminal_ineq_constraint_a1;
g{T+1,2} = @terminal_ineq_constraint_b1;
g{T+1,3} = @terminal_ineq_constraint_c1;
g{T+1,4} = @terminal_ineq_constraint_c2;

h{T+1,1} = @terminal_constraint_a1;
h{T+1,2} = @terminal_constraint_b1;
h{T+1,3} = @terminal_constraint_c1;
h{T+1,4} = @terminal_constraint_c2;

%% Options
params.tauval_init = {.1,.1,.1,.1,.1,.1};
params.tauval_decrease = 0.1;
params.tauval_tolerance = 0.001;
params.resid_tolerance = 0.001;
params.debug_plot = true;
params.debug_print = 2;
params.max_linesearch_iters = 3;
params.warm_up_iters = 0;
params.warm_up_alpha = 0.01;
params.max_cl_linesearch_iters = 5;
params.alpha_init = 1;
params.failed_ls_alpha = 1;
params.beta = 0.5;
params.open_loop = false;
params.wait_for_ol_convergence = true;
params.crash_start = true;
params.independent_updates = true;
params.feasible = false;
params.use_initialization = false;
params.max_feasible_iters = 30;
params.constraint_tolerance = 1e-3;
params.tauval_max = 1;

%% Run 
prob_c1 = 1;
[residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);
solution_set{1} = xval;

params.use_initialization = true;
params.init_step = 0;
params.wait_for_ol_convergence = false;
params.crash_start = false;
params.tauval_init = {.01,.01,.01,.01,.01,.01};

prob_c1 = .667;
[residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
solution_set{2} = xval;
prob_c1 = .333;
[residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
solution_set{3} = xval;
prob_c1 = 0;
[residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
solution_set{4} = xval;
plot_solution(xval,x0);
%%
% params.tauval_init = {.01,.1,.01,.1,.01,.1};
% params.crash_start = false;
% params.wait_for_ol_convergence = true;
% params.use_initialization = true;
% params.warm_up_iters = 0;
% params.debug_plot = true;
% 
% 
% for iter = 1:20
%     mpc_iter = iter
%     x0 = xval{20};
%     params.init_step = 20;
%     x0(5:8) = x0(1:4);
%     x0(13:16) = x0(9:12);
%     x0(21:24) = x0(17:20);
%     [residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution, initialization);
%     solution_set{iter+1} = xval;
% %     plot_solution(xval,x0);
% end
%%
disp('done');
% plot_solution(solution_set{8},solution_set{8}{1});
plot_solution_set(solution_set);

%% Plotting

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

    for t = 1:T+1
        xx(t,:) = xval{t};
        if plot_polite
            xxp(t,:) = varargin{1}{t};
        end
        
        spot_a1 = viscircles(xx(t,1:2),avoid_r/2,'Color', 'b');
        car_a1 = draw_car(xx(t,1:2),'b');

        spot_b1 = viscircles(xx(t,5:6),avoid_r/2,'Color', 'r');
        car_b1 = draw_car(xx(t,5:6),'r');

        spot_c1 = viscircles(xx(t,9:10),avoid_r/2,'Color', 'g');
        car_c1 = draw_car(xx(t,9:10),'g');
        spot_c2 = viscircles(xx(t,13:14),avoid_r/2,'Color', 'c');
        car_c2 = draw_car(xx(t,13:14),'c');

%         if plot_polite
%             pspot_a = viscircles([xxp(t+1,1),xxp(t+1,2)],avoid_r/2,'Color', 'b','LineStyle',':');
%             pcar_a = patch('XData',[xxp(t+1,1)-half_w,xxp(t+1,1)-half_w,xxp(t+1,1)+half_w,xxp(t+1,1)+half_w,xxp(t+1,1)-half_w],...
%                    'YData',[xxp(t+1,2)-half_l,xxp(t+1,2)+half_l,xxp(t+1,2)+half_l,xxp(t+1,2)-half_l,xxp(t+1,2)-half_l],...
%                        'FaceAlpha',1,'FaceColor', 'b');
%             pspot_b = viscircles([xxp(t+1,5),xxp(t+1,6)],avoid_r/2,'Color', 'r','LineStyle',':');
%             pcar_b = patch('XData',[xxp(t+1,5)-half_w,xxp(t+1,5)-half_w,xxp(t+1,5)+half_w,xxp(t+1,5)+half_w,xxp(t+1,5)-half_w],...
%                    'YData',[xxp(t+1,6)-half_l,xxp(t+1,6)+half_l,xxp(t+1,6)+half_l,xxp(t+1,6)-half_l,xxp(t+1,6)-half_l],...
%                        'FaceAlpha',1,'FaceColor', 'r');
%             pspot_c1 = viscircles([xxp(t+1,9),xxp(t+1,10)],avoid_r/2,'Color', 'g','LineStyle',':');
%             pcar_c1 = patch('XData',[xxp(t+1,9)-half_w,xxp(t+1,9)-half_w,xxp(t+1,9)+half_w,xxp(t+1,9)+half_w,xxp(t+1,9)-half_w],...
%                        'YData',[xxp(t+1,10)-half_l,xxp(t+1,10)+half_l,xxp(t+1,10)+half_l,xxp(t+1,10)-half_l,xxp(t+1,10)-half_l],...
%                        'FaceAlpha',1,'FaceColor', 'g');
%             pspot_c2 = viscircles([xxp(t+1,13),xxp(t+1,14)],avoid_r/2,'Color', 'c','LineStyle',':');
%             pcar_c2 = patch('XData',[xxp(t+1,13)-half_w,xxp(t+1,13)-half_w,xxp(t+1,13)+half_w,xxp(t+1,13)+half_w,xxp(t+1,13)-half_w],...
%                        'YData',[xxp(t+1,14)-half_l,xxp(t+1,14)+half_l,xxp(t+1,14)+half_l,xxp(t+1,14)-half_l,xxp(t+1,14)-half_l],...
%                        'FaceAlpha',1,'FaceColor', 'g');
% 
%         end
        
        axis([xx(t,9)-20,xx(t,9)+20,xx(t,10)-20,xx(t,10)+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.01);
        end
        delete(spot_a1);
        delete(car_a1);
        delete(spot_b1);
        delete(car_b1);
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
    if make_vid
        close(writerObj);
    end
end

%%
function plot_solution_set(solution_set)
    make_vid = true;
    plot_polite = false;
    close all;
    fh = figure; 

    if make_vid
        writerObj = VideoWriter('polite_merge','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    hold on;
    global border_left border_right avoid_r linear_dyn T half_l half_w;
    
    plot([border_left, border_left],[-200,200],'-k');
    plot([border_left/2, border_left/2],[-200,200],'--k');
    plot([0, 0],[-200,200],'-k');
    plot([border_right/2, border_right/2],[-200,200],'--k');
    plot([border_right, border_right],[-200,200],'-k');
    axis([-20,20,-20,20]);

    
    L = length(solution_set);
    green = [0,1,0];
    yellow = [1,1,0];
    blue = [0,0,1];
    probs = [0,.333,.667,1];
    for t = 1:T+1
        for s = 1:L
            xx(t,:) = solution_set{s}{t};
            color = 'b';
            spot_a1{s} = viscircles(xx(t,1:2),avoid_r/2,'Color', color);
            car_a1{s} = draw_car(xx(t,1:2),color);
    %             spot_a2 = viscircles(xx(ind,5:6),avoid_r/2,'Color', 'b','LineStyle',':');
    %             car_a2 = draw_car(xx(ind,5:6),'b');

            spot_b1{s} = viscircles(xx(t,5:6),avoid_r/2,'Color', 'r');
            car_b1{s} = draw_car(xx(t,5:6),'r');
    %             spot_b2 = viscircles(xx(ind,13:14),avoid_r/2,'Color', 'r','LineStyle',':');
    %             car_b2 = draw_car(xx(ind,13:14),'r');
            color = 'g';
            spot_c1{s} = viscircles(xx(t,9:10),avoid_r/2,'Color', color);
            car_c1{s} = draw_car(xx(t,9:10),color);
            spot_c2{s} = viscircles(xx(t,13:14),avoid_r/2,'Color', color);
            car_c2{s} = draw_car(xx(t,13:14),color);
        end

        axis([-20,20,xx(t,2)-20,xx(t,2)+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.01);
        end
        for s = 1:L
            delete(spot_a1{s});
            delete(car_a1{s});

            delete(spot_b1{s});
            delete(car_b1{s});

            delete(spot_c1{s});
            delete(car_c1{s});
            delete(spot_c2{s});
            delete(car_c2{s});
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
    global avoid_r;
    
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    
    const = [(xa1-xb1)'*(xa1-xb1)-avoid_r*avoid_r;
             right_lane_ineq_constraint(xa1);
             left_lane_ineq_constraint(xa1)];
end

function const = running_ineq_constraint_b1(x)
    global avoid_r;
    
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
   
    const = [right_lane_ineq_constraint(xb1);
             left_lane_ineq_constraint(xb1)];
end

function const = running_ineq_constraint_c1(x)
    global avoid_r;
    
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    const = [(xc1-xb1)'*(xc1-xb1)-avoid_r*avoid_r;
             (xc1-xa1)'*(xc1-xa1)-avoid_r*avoid_r];
    const = [const;
            right_lane_ineq_constraint(xc1);
            left_lane_ineq_constraint(xc1)];
end

function const = running_ineq_constraint_c2(x)
    global avoid_r;
    
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    xc2 = x(13:14);
    const = [(xc2-xb1)'*(xc2-xb1)-avoid_r*avoid_r;
             (xc2-xa1)'*(xc2-xa1)-avoid_r*avoid_r];
    const = [const;
             right_lane_ineq_constraint(xc2);
             left_lane_ineq_constraint(xc2)];
end

%% Terminal ineq Constraints

function const = terminal_ineq_constraint_a1(x)
    const = [running_ineq_constraint_a1(x)];
end

function const = terminal_ineq_constraint_b1(x)    
    const = [running_ineq_constraint_b1(x)];
end

function const = terminal_ineq_constraint_c1(x)
    const = [running_ineq_constraint_c1(x)];
end

function const = terminal_ineq_constraint_c2(x)
    const = [running_ineq_constraint_c2(x)];
end

%% Terminal constraints 

function const = terminal_constraint_a1(x)
    xa1 = x(1:2);
    const = left_lane_constraint(xa1);
end


function const = terminal_constraint_b1(x)
    xb1 = x(5:6);
    const = right_lane_constraint(xb1);
end

function const = terminal_constraint_c1(x)
    xc1 = x(9:10);
    const = left_lane_constraint(xc1);
end

function const = terminal_constraint_c2(x)
    xc2 = x(13:14);
    const = left_lane_constraint(xc2);
end

function val = terminal_cost(~)
    val = 0;
end

%% Running costs
function val = running_cost_a1(xu)
    global gv_a1;
    global center_lane_pen speed_pen;
    xa1 = xu(1:4);
    ua1 = xu(17:18);
    val = ua1'*ua1 + speed_pen*speed_cost(xa1,gv_a1)+center_lane_pen*center_lane_cost(xa1);
end

function val = running_cost_b1(xu)
    global gv_b1;
    global center_lane_pen speed_pen;
    xb1 = xu(5:9);
    ub1 = xu(19:20);
    val = ub1'*ub1 + speed_pen*speed_cost(xb1,gv_b1)+center_lane_pen*center_lane_cost(xb1);
end

function val = running_cost_c1(xu)
    global gv_c1;
    global center_lane_pen speed_pen;
    xc1 = xu(9:12);
    uc1 = xu(21:22);
    val = uc1'*uc1 + speed_pen*speed_cost(xc1,gv_c1)+center_lane_pen*center_lane_cost(xc1);
end

function val = running_cost_c2(xu)
    global gv_c2;
    global center_lane_pen speed_pen;
    xc2 = xu(13:16);
    uc2 = xu(23:24);
    val = uc2'*uc2 + speed_pen*speed_cost(xc2,gv_c2)+center_lane_pen*center_lane_cost(xc2);
end

%% Polite costs
function val = running_cost_a1_polite(xu)
    global prob_c1;
    val = running_cost_a1(xu)+prob_c1*running_cost_c1(xu)+(1-prob_c1)*running_cost_c2(xu);
end

function val = running_cost_b1_polite(xu)
    val = running_cost_b1(xu);
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

function val = center_lane_cost(x)
    global border_left border_right;
    val = ((x(1)-border_left/2)*(x(1)-border_right/2))^2;
end

function const = left_lane_constraint(x)
    global border_left;
    const = x(1)-border_left/2;
end

function const = right_lane_constraint(x)
    global border_right;
    const = x(1)-border_right/2;
end

function const = left_lane_ineq_constraint(x)
    global border_left half_w;
    const = x(1)-half_w-border_left;
end

function const = right_lane_ineq_constraint(x)
    global border_right half_w;
    const = border_right-half_w-x(1);
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
    xb1 = xu(5:8);
    xc1 = xu(9:12);
    xc2 = xu(13:16);
    
    ua1 = xu(17:18);
    ub1 = xu(19:20);
    uc1 = xu(21:22);
    uc2 = xu(23:24);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    zc1 = xc1 + dt*unicycle_dyn(xc1,uc1);
    zc2 = xc2 + dt*unicycle_dyn(xc2,uc2);
    next = [za1;zb1;zc1;zc2];
end