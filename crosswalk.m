%% crosswalk game
clear; clc; close all;

global dt;
dt = 0.2;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;

border_left = -4;
border_right = 4;

global center_lane_pen speed_pen lane_pen;
center_lane_pen = 0;
speed_pen = 1;
lane_pen = 5;

global car_r ped_r avoid_r half_l half_w;
half_l = 1.5;
half_w = 0.75;
car_r = sqrt(half_l^2+half_w^2);
ped_r = 2;
avoid_r = (ped_r+car_r);

global ped_w;
ped_w = 0.5;

global crosswalk_position crosswalk_length num_stripes;
crosswalk_position = 0;
crosswalk_length = 2;
num_stripes = 6;

global prob_b1;

prob_b1 = 0;

global gv_a1 gv_b1 gv_b2;



global T
T = 30;
needed_speed = -10/(T*dt);


gv_a1 = -.5;
gv_b1 = 1;
gv_b2 = 2;

x0 = [border_right+1; crosswalk_position+crosswalk_length/2;gv_a1;0;
      border_right/2;-18;0;gv_b1;
      border_right/2;-18;0;gv_b2];
  

rollout_x0 = [border_right+1; crosswalk_position+crosswalk_length/2;gv_a1;0;
      border_right/2;-18;0;gv_b1;
      border_right/2;-18;0;gv_b1];
  
global unimpeded;

unimpeded = x0(6)+T*dt*gv_b2;

m{1} = 2;
m{2} = 2;
m{3} = 2;
N = 3;
n = N*4;


for t = 1:T
    l{t,1} = @running_cost_a1_polite;
    l{t,2} = @running_cost_b1;
    l{t,3} = @running_cost_b2;
  
    g{t,1} = @running_ineq_constraint_a1;
    g{t,2} = @running_ineq_constraint_b1;
    g{t,3} = @running_ineq_constraint_b2;
    
%     g{t,1} = @empty_constraint;
%     g{t,2} = @empty_constraint;
%     g{t,3} = @empty_constraint;

    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;


end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
l{T+1,3} = @empty_cost;

g{T+1,1} = @terminal_ineq_constraint_a1;
g{T+1,2} = @terminal_ineq_constraint_b1;
g{T+1,3} = @terminal_ineq_constraint_b2;

% g{T+1,1} = @empty_constraint;
% g{T+1,2} = @empty_constraint;
% g{T+1,3} = @empty_constraint;

h{T+1,1} = @empty_constraint;
h{T+1,2} = @empty_constraint;
h{T+1,3} = @empty_constraint;

%% Options
params.tauval_init = {1,1,1};
params.tauval_decrease = 0.1;
params.tauval_tolerance = 0.003;
params.resid_tolerance = 0.01;
params.debug_plot = true;
params.debug_print = 2;
params.max_linesearch_iters = 10;
params.warm_up_iters = 0;
params.warm_up_alpha = 0.01;
params.max_cl_linesearch_iters = 10;
params.alpha_init = 1;
params.failed_ls_alpha = 1;
params.beta = 0.5;
params.open_loop = false;
params.crash_start = false;
params.independent_updates = false;
params.feasible = false;
params.use_initialization = false;
params.max_feasible_iters = 100;
params.constraint_tolerance = 1e-3;
params.tauval_max = 1;
params.rollout_x0 = rollout_x0;
params.use_rollout_x0 = true;
params.converge_on_every_policy = true;
params.wait_for_ol_convergence = true;
params.basic_tauschedule = true;
params.fresh_policy_threshold = 3;
params.slack_min = 1e-4;
params.restore_after_failed_ls = false;
params.dev_penalty = 30;
params.dev_decay = 0.5;
params.single_alpha = true;
params.frac_to_bound = 0.995;


%% Run 
prob_b1 = 1;
params.open_loop = false;
[residuals, xval, initialization] = ec_solver_clean(@f, h, g, l, n, m, N, T, x0, params, @plot_solution, 0);
solution_set{1} = xval;
% 
% params.use_initialization = true;
% params.init_step = 0;
% params.open_loop = false;
% params.wait_for_ol_convergence = true;
% params.crash_start = false;
% params.tauval_init = {.1,.1,.1};
% params.independent_updates = false;
% 
prob_b1 = .1
[residuals, xval, initialization2] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
solution_set{2} = xval;
%%
plot_solution_set(solution_set)

% prob_b1 = .33
% [residuals, xval, initialization2] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
% solution_set{3} = xval;
% 
% 
% prob_b1 = 0
% [residuals, xval, initialization2] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
% solution_set{4} = xval;
% prob_b1 = 0;
% [residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution,initialization);
% solution_set{4} = xval;
% plot_solution(xval,x0);
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
% disp('done');
% plot_solution(solution_set{8},solution_set{8}{1});
% plot_solution_set(solution_set);

%% Plotting

function plot_solution(xval,x0)
    make_vid = false;
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
    global car_r ped_r T unimpeded;
    
    draw_lanes();
    draw_crosswalk();
    plot([-10,10],[unimpeded,unimpeded],'--r');

    for t = 1:T+1
        xx(t,:) = xval{t};
        
        spot_a1 = viscircles(xx(t,1:2),ped_r,'Color', 'b','LineStyle',':');
        ped_a1 = draw_ped(xx(t,1:2),'b');

        spot_b1 = viscircles(xx(t,5:6),car_r,'Color', 'g','LineStyle',':');
        car_b1 = draw_car(xx(t,5:6),'r');
        
        spot_b2 = viscircles(xx(t,9:10),car_r,'Color', 'g','LineStyle',':');
        car_b2 = draw_car(xx(t,9:10),'r');
        
        axis([-20,+20,xx(t,6)-20,xx(t,6)+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(.01);
        end
        delete(spot_a1);
        delete(ped_a1);
        delete(spot_b1);
        delete(car_b1);
        delete(spot_b2);
        delete(car_b2);
        
    end
    if make_vid
        close(writerObj);
    end
end

%%
function plot_solution_set(solution_set)
    make_vid = false;
    plot_polite = false;
    close all;
    fh = figure; 

    if make_vid
        writerObj = VideoWriter('kiwibot','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
%         axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    hold on;
    global car_r ped_r T;
    
    draw_lanes();
    draw_crosswalk();
    
    L = length(solution_set);
    for t = 1:T+1
        for s = 1:L
            xx(t,:) = solution_set{s}{t};
            spot_a1{s} = viscircles(xx(t,1:2),ped_r,'Color', 'b','LineStyle',':');
            ped_a1{s} = draw_ped(xx(t,1:2),'b');

            spot_b1{s} = viscircles(xx(t,5:6),car_r,'Color', 'g','LineStyle',':');
            car_b1{s} = draw_car(xx(t,5:6),'r');

            spot_b2{s} = viscircles(xx(t,9:10),car_r,'Color', 'g','LineStyle',':');
            car_b2{s} = draw_car(xx(t,9:10),'r');

            
        end
        axis([-20,+20,-20,+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.01);
        end
        for s = 1:L
            delete(spot_a1{s});
            delete(ped_a1{s});
            delete(spot_b1{s});
            delete(car_b1{s});
            delete(spot_b2{s});
            delete(car_b2{s});
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
    xb2 = x(9:10);
    const = zeros(0,1);
%     const = [crosswalk_constraint(xa1)];
end

function const = running_ineq_constraint_b1(x)
    global avoid_r;
    
    xa1 = x(1:2);
    xb1 = x(5:6);
    xb2 = x(9:10);
  
    const = [(xb1-xa1)'*(xb1-xa1)-avoid_r*avoid_r];
%              right_lane_ineq_constraint(xb1);
%              left_lane_ineq_constraint(xb1)];
end

function const = running_ineq_constraint_b2(x)
    global avoid_r;
    
    xa1 = x(1:2);
    xb1 = x(5:6);
    xb2 = x(9:10);

    const = [(xb2-xa1)'*(xb2-xa1)-avoid_r*avoid_r];
%             right_lane_ineq_constraint(xb2);
%             left_lane_ineq_constraint(xb2)];
end



%% Terminal ineq Constraints

function const = terminal_ineq_constraint_a1(x)
    const = [running_ineq_constraint_a1(x);
             terminal_constraint_a1(x)];
end

function const = terminal_ineq_constraint_b1(x)    
    const = [running_ineq_constraint_b1(x)];
end

function const = terminal_ineq_constraint_b2(x)
    const = [running_ineq_constraint_b2(x)];
end

%% Terminal constraints 

function const = terminal_constraint_a1(x)
    xa1 = x(1:2);
    const = crosswalk_cross_left_constraint(xa1);
end


function const = terminal_constraint_b1(x)
    xb1 = x(5:6);
    const = right_lane_constraint(xb1);
end

function const = terminal_constraint_b2(x)
    xb2 = x(9:10);
    const = right_lane_constraint(xb2);
end

function val = terminal_cost(~)
    val = 0;
end

%% Running costs
function val = running_cost_a1(xu)
    global gv_a1 speed_pen;
    xa1 = xu(1:4);
    ua1 = xu(13:14);
    val = ua1'*ua1 + speed_pen*(xa1(3)-gv_a1)^2+crosswalk_cost(xa1);
end

function val = running_cost_b1(xu)
    global gv_b1;
    global speed_pen lane_pen;
    
    R = diag([10,1]);
    xb1 = xu(5:8);
    ub1 = xu(15:16);
    val = ub1'*R*ub1 + speed_pen*speed_cost(xb1,gv_b1) + lane_pen*right_lane_cost(xb1);
end

function val = running_cost_b2(xu)
    global gv_b2;
    global speed_pen lane_pen;
    R = diag([10,1]);
    xb2 = xu(9:12);
    ub2 = xu(17:18);
    val = ub2'*R*ub2 + speed_pen*speed_cost(xb2,gv_b2) + lane_pen*right_lane_cost(xb2);
end


%% Polite costs
function val = running_cost_a1_polite(xu)
    global prob_b1;
    val = running_cost_a1(xu)+prob_b1*running_cost_b1(xu)+(1-prob_b1)*running_cost_b2(xu);
end

function val = running_cost_b1_polite(xu)
    val = running_cost_b1(xu);
end

function val = running_cost_b2_polite(xu)
    val = running_cost_b2(xu);
end
%% Helpers

function val = speed_cost(x,goal)
    val = (x(4)-goal)^2;
end

function val = either_lane_cost(x)
    global border_left border_right;
    val = ((x(1)-border_left/2)*(x(1)-border_right/2))^2;
end

function val = right_lane_cost(x)
    global border_right;
    val = (x(1)-border_right/2)^2;
end

function val = left_lane_cost(x)
    global border_left;
    val = (x(1)-border_left/2)^2;
end

function val = crosswalk_cost(x)
    global crosswalk_length crosswalk_position;
    val = (x(2)-crosswalk_position-crosswalk_length/2)^2;
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

function const = crosswalk_constraint(x)
    global crosswalk_length crosswalk_position ped_w;
    const = [crosswalk_position+crosswalk_length+5-(x(2)+ped_w);
             x(2)-ped_w - crosswalk_position+5];
end

function const = crosswalk_cross_left_constraint(x)
    global border_left;
    const = [border_left-1-x(1)];
end

function pt = draw_car(center,color)
    global half_l half_w;
    pt = patch('XData',[center(1)-half_w,center(1)-half_w,center(1)+half_w,center(1)+half_w,center(1)-half_w],...
               'YData',[center(2)-half_l,center(2)+half_l,center(2)+half_l,center(2)-half_l,center(2)-half_l],...
               'FaceAlpha',1,'FaceColor', color);
end

function pt = draw_ped(center,color)
    global ped_w;
    pt = patch('XData',[center(1)-ped_w,center(1)-ped_w,center(1)+ped_w,center(1)+ped_w,center(1)-ped_w],...
               'YData',[center(2)-ped_w,center(2)+ped_w,center(2)+ped_w,center(2)-ped_w,center(2)-ped_w],...
               'FaceAlpha',1,'FaceColor', color);
end

function draw_lanes()
    global border_left border_right;
    plot([border_left, border_left],[-200,200],'-k');
    plot([border_left/2, border_left/2],[-200,200],'--k');
    plot([0, 0],[-200,200],'-k');
    plot([border_right/2, border_right/2],[-200,200],'--k');
    plot([border_right, border_right],[-200,200],'-k');
end

function draw_crosswalk()
    global crosswalk_position crosswalk_length border_left border_right num_stripes;
    stripe_width = (border_right-border_left) / (2*num_stripes+1);
    loc = border_left+stripe_width;
    for i = 1:num_stripes
        draw_stripe([loc,crosswalk_position],stripe_width, crosswalk_length);
        loc = loc + 2*stripe_width;
    end
end

function draw_stripe(location, width, length)
    x = location(1);
    y = location(2);
    pt = patch('XData',[x,x,x+width,x+width,x],...
               'YData',[y,y+length,y+length,y,y],...
               'FaceAlpha',1,'FaceColor', 'w');
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
    xb2 = xu(9:12);
    
    ua1 = xu(13:14);
    ub1 = xu(15:16);
    ub2 = xu(17:18);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    zb2 = xb2 + dt*unicycle_dyn(xb2,ub2);
    next = [za1;zb1;zb2];
end