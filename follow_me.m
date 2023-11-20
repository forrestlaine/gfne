%% follow-me nonlinear dynamics game
clear; clc; close all;

global dt;
dt = 0.1;
global T
T = 100;

global linear_dyn;
linear_dyn = true;

global border_left border_right border_top border_bottom;


border_left = -10;
border_right = 10;
border_top = 10;
border_bottom = -10;

global QQ;

qq = [1/(border_right-border_left); 1/(border_top-border_bottom); T*dt/(border_top-border_bottom); 1/(2*pi)];
QQ = diag(qq);

global RR;
rr = [1/(border_right-border_left); 1/(2*pi)];
RR = diag(rr);


global car_r ped_r avoid_r half_l half_w;
half_l = 1.5;
half_w = 0.75;
car_r = sqrt(half_l^2+half_w^2);
ped_r = 2;
avoid_r = (car_r+car_r);





needed_speed = -10/(T*dt);

global xga xgb;

x0 = [-9;-3;1;pi/3;
      -6; 6;1;-pi/6];
  
rollout_x0 = x0;
% rollout_x0(3) = 0.1;
% rollout_x0(7) = 1.5;
  
  
xga = [6;6;1;pi/3];
xgb = [6;-6;1;-pi/3];
  


m{1} = 2;
m{2} = 2;
N = 2;
n = N*4;


for t = 1:T
    l{t,1} = @running_cost_a1;
    l{t,2} = @running_cost_b1_polite;
  
    g{t,1} = @running_ineq_constraint_a1;
    g{t,2} = @running_ineq_constraint_b1;
    
%     g{t,1} = @empty_constraint;
%     g{t,2} = @empty_constraint;
%     g{t,3} = @empty_constraint;

    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;


end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;
% l{T+1,1} = @terminal_cost_a1;
% l{T+1,2} = @terminal_cost_b1;

g{T+1,1} = @terminal_ineq_constraint_a1;
g{T+1,2} = @terminal_ineq_constraint_b1;

g{T+1,1} = @empty_constraint;
g{T+1,2} = @empty_constraint;
% g{T+1,3} = @empty_constraint;

h{T+1,1} = @terminal_constraint_a1;
h{T+1,2} = @terminal_constraint_b1;
% h{T+1,1} = @empty_constraint;
% h{T+1,2} = @empty_constraint;
%% Options
params.tauval_init = {.1,.1,.1};
params.tauval_decrease = 0.1;
params.tauval_tolerance = 0.001;
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
params.converge_on_every_policy = false;
params.wait_for_ol_convergence = false;
params.basic_tauschedule = false;
params.fresh_policy_threshold = 2;
params.slack_min = 1e-4;
params.restore_after_failed_ls = false;
params.dev_penalty = 1;
params.dev_decay = 0.5;
params.single_alpha = true;


%% Run 
params.open_loop = false;
[residuals, xval, initialization] = ec_solver_ip(@f, h, g, l, n, m, N, T, x0, params, @plot_solution);
solution_set{1} = xval;


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
    global car_r T;
    
    draw_borders();

    for t = 1:T+1
        xx(t,:) = xval{t};
        
        spot_a1 = viscircles(xx(t,1:2),car_r,'Color', 'b','LineStyle',':');
        car_a1 = draw_car(xx(t,1:2),xx(t,4),'b');

        spot_b1 = viscircles(xx(t,5:6),car_r,'Color', 'g','LineStyle',':');
        car_b1 = draw_car(xx(t,5:6),xx(t,8),'r');
        
        axis([-10,+10,-10,+10]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(.01);
        end
        delete(spot_a1);
        delete(car_a1);
        delete(spot_b1);
        delete(car_b1);
        
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
    xa1 = x(1:2);
    xb1 = x(5:6);
    const = avoid_constraint(xa1,xb1);
end

function const = running_ineq_constraint_b1(x)
    
    xa1 = x(1:2);
    xb1 = x(5:6);
%     const = avoid_constraint(xa1,xb1);
    const = zeros(0,1);
end



%% Terminal ineq Constraints

function const = terminal_ineq_constraint_a1(x)
    const = [running_ineq_constraint_a1(x)];
end

function const = terminal_ineq_constraint_b1(x)    
    const = [running_ineq_constraint_b1(x)];
end


%% Terminal constraints 

function const = terminal_constraint_a1(x)
    global xga;
    xa1 = x(1:4);
    xb1 = x(5:8);
    diff = xga-xa1;
    const = diff;
end


function const = terminal_constraint_b1(x)
    global xgb;
    xa1 = x(1:4);
    xb1 = x(5:8);
    diff = xgb-xb1;
    const = diff;
end

function val = terminal_cost_a1(x)
    global QQ;
    xa1 = x(1:4);
    xb1 = x(5:8);
    val = 0;
    
    val = val+100*xb1(1)*QQ(1,1)*xb1(1);
    val = val+100*xb1(2)*QQ(2,2)*xb1(2);
    val = val+100*xb1(3)*QQ(3,3)*xb1(3);
    angle = xb1(4);
    angle = mod(angle+pi,2*pi)-pi;
    val = val+100*angle*QQ(4,4)*angle;
end

function val = terminal_cost_b1(x)
    global QQ;
    xa1 = x(1:4);
    xb1 = x(5:8);
    val = 0;
    xd = xa1-xb1;
    val = val+100*xd(1)*QQ(1,1)*xd(1);
    val = val+100*xd(2)*QQ(2,2)*xd(2);
    val = val+100*xd(3)*QQ(3,3)*xd(3);
    angle = xd(4);
    angle = mod(angle+pi,2*pi)-pi;
    val = val+100*angle*QQ(4,4)*angle;
end

%% Running costs
function val = running_cost_a1(xu)
    global RR;
    xa1 = xu(1:4);
    ua1 = xu(9:10);
    val = ua1'*RR*ua1;
end

function val = running_cost_b1(xu)
    global RR;
    xb1 = xu(5:8);
    ub1 = xu(11:12);
    val = ub1'*RR*ub1;
end

%% Polite costs
function val = running_cost_a1_polite(xu)
    val = running_cost_a1(xu)+running_cost_b1(xu);
end

function val = running_cost_b1_polite(xu)
    val = running_cost_b1(xu)+running_cost_a1(xu);
end

%% Helpers


function const = boundary_constraint(x)
    global border_left border_right border_top border_bottom car_r;
    const = [border_right - (x(1)+car_r);
             (x(1)-car_r) - border_left;
             border_top - (x(2)+car_r);
             (x(2)-car_r) - border_bottom];
end

function const = avoid_constraint(x1,x2)
    global avoid_r;
    const = (x1-x2)'*(x1-x2)-avoid_r^2;
end

function pt = draw_car(center,angle,color)
    global half_l half_w;
    pt = patch('XData',[center(1)-half_w,center(1)-half_w,center(1)+half_w,center(1)+half_w,center(1)-half_w],...
               'YData',[center(2)-half_l,center(2)+half_l,center(2)+half_l,center(2)-half_l,center(2)-half_l],...
               'FaceAlpha',1,'FaceColor', color);
    rotate(pt, [0,0,1],90+360/(2*pi)*angle,[center,0]);
end

function draw_borders()
    global border_left border_right border_top border_bottom;
    rectangle('Position',[border_left,border_bottom,border_right-border_left,border_top-border_bottom]);
end

function vec = unicycle_dyn(x,u)
    vec = x(3)*cos(x(4));
    vec = [vec; x(3)*sin(x(4))];
    vec = [vec; u(1)];
    vec = [vec; u(2)];
end

function next = f(xu)
    global dt;
    xa1 = xu(1:4);
    xb1 = xu(5:8);
    
    ua1 = xu(9:10);
    ub1 = xu(11:12);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    next = [za1;zb1];
end