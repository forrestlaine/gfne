%% Is there a car?
clear; clc; close all;

global dt;
dt = 0.5;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;

border_left = -4;
border_right = 4;

global center_lane_pen speed_pen lane_pen;
center_lane_pen = 0;
speed_pen = 100;
lane_pen = 5;
global R;
R = diag([.1,.1]);

global car_r ped_r avoid_r half_l half_w;
half_l = 1.5;
half_w = 0.75;
car_r = sqrt(half_l^2+half_w^2);
ped_r = 2;
avoid_r = (car_r+car_r);

global ped_w;
ped_w = 0.5;

global crosswalk_position crosswalk_length num_stripes;
crosswalk_position = 0;
crosswalk_length = 2;
num_stripes = 6;

global prob_b1 prob_b2;


global gv_a1 gv_b1 gv_b2 gv_c1;



global T
T = 50;
needed_speed = -10/(T*dt);


gv_a1 = 1;
gv_b1 = 2.0;
gv_b2 = 3.0;
gv_c1 = 1;

x0 = [border_right/2; 0;0; gv_a1;
      border_left/2; -10;0;gv_b1;
      border_left/2; -10;0;gv_b2;
      border_right/2; -5;0;gv_c1];
  

rollout_x0 = [border_right/2; 0; -0.1; gv_a1;
      border_left/2;  -10;0;0;
      border_left/2; -10;0;0;
      border_right/2; -5;0;gv_c1];
  
m{1} = 2;
m{2} = 2;
m{3} = 2;
m{4} = 2;
N = 4;
n = N*4;


for t = 1:T
    if t > 20
        l{t,1} = @running_cost_a1_polite_ll;
    else
        l{t,1} = @running_cost_a1_polite;
    end
    l{t,2} = @running_cost_b1;
    l{t,3} = @running_cost_b2;
    l{t,4} = @running_cost_c1;
  
    g{t,1} = @running_ineq_constraint_a1;
    g{t,2} = @running_ineq_constraint_b1;
    g{t,3} = @running_ineq_constraint_b2;
    g{t,4} = @running_ineq_constraint_c1;
    


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
g{T+1,3} = @terminal_ineq_constraint_b2;
g{T+1,4} = @terminal_ineq_constraint_c1;

h{T+1,1} = @terminal_constraint_a1;
h{T+1,2} = @terminal_constraint_b1;
h{T+1,3} = @terminal_constraint_b2;
h{T+1,4} = @terminal_constraint_c1;


%% Run 
prob_b1 = .99;
prob_b2 = .01;
params.use_ws = false;
params.open_loop = true;
params.contact = false;
params.use_initialization = false;
params.debug_plot = false;
params.linear = true;
[xval,ws,durs] = fb_game_solver_active_set(@f, h, g, l, n, m, N, T, x0, rollout_x0, params, @plot_solution);
solution_set{1} = xval;
durs.p = prob_b1;
durations{1} = durs;

prob_b2 = 2/3;
prob_b1 = 1/3;
params.use_ws = false;
params.use_initialization = false;
[xval,ws,durs] = fb_game_solver_active_set(@f, h, g, l, n, m, N, T, x0, rollout_x0, params, @plot_solution,ws,xval);
solution_set{2} = xval;
durs.p = prob_b1;
durations{2} = durs;

prob_b2 = .99;
prob_b1 = 0.01;
[xval,ws,durs] = fb_game_solver_active_set(@f, h, g, l, n, m, N, T, x0, rollout_x0, params, @plot_solution,ws,xval);
solution_set{3} = xval;
durs.p = prob_b1;
durations{3} = durs;


%%
% for i = 1:5
%     avg_speed = 0;
%     max_speed = 0;
%     for t=1:101
%         spd = solution_set{i}{t}(4);
%         avg_speed = avg_speed+spd;
%         max_speed = max(max_speed,spd);
%     end
%     average_speeds(i) = avg_speed/101;
%     max_speeds(i) = max_speed;
% end
%         
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
% plot_solution_set({solution_set{1},solution_set{2},solution_set{3}});

%%


%% Plotting

function plot_solution(xval,x0)
    make_vid = false;
    close all;
    fh = figure; 
    hold on;
    if make_vid
        writerObj = VideoWriter('drive_slow','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    global car_r T;
    
    draw_lanes();

    for t = 1:T+1
        xx(t,:) = xval{t};
        
       
        spot_b1 = viscircles(xx(t,5:6),car_r,'Color', [1,.7,0],'LineStyle',':');
        car_b1 = draw_car(xx(t,5:6),[1,.7,0],1);
        
        spot_b2 = viscircles(xx(t,9:10),car_r,'Color', 'r','LineStyle',':');
        car_b2 = draw_car(xx(t,9:10),'r',1);
        
        spot_a1 = viscircles(xx(t,1:2),car_r,'Color', 'b','LineStyle',':');
        car_a1 = draw_car(xx(t,1:2),'b',1);
        
        spot_c1 = viscircles(xx(t,13:14),car_r,'Color', 'g','LineStyle',':');
        car_c1 = draw_car(xx(t,13:14),'g',1);
        
        axis([-20,+20,xx(t,6)-20,xx(t,6)+20]);
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
        delete(spot_b2);
        delete(car_b2);
        delete(spot_c1);
        delete(car_c1);
        
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
        writerObj = VideoWriter('merge_in_front_spedup','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
%         axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    hold on;
    global car_r ped_r T;
    
    draw_lanes();
%     draw_crosswalk();
    alphas = [.2,.6,1];
    
    L = length(solution_set);
    for t = 1:T+1
%             spot_b1{s} = viscircles(xx(t,5:6),car_r,'Color', alphas(s)*[1,0,0],'LineStyle',':');
%         for s = 1:L
%             xx(t,:) = solution_set{s}{t};    
%             car_b1{s} = draw_car(xx(t,5:6),[1,.7,0],alphas(s));
%         end
        for s = 1:L
            xx(t,:) = solution_set{s}{t};
%             spot_b2{s} = viscircles(xx(t,9:10),car_r,'Color', 'r','LineStyle',':');
            car_b2{s} = draw_car(xx(t,9:10),[1,0,0],alphas(s));
        end
        for s = 1:L
            xx(t,:) = solution_set{s}{t};
            car_a1{s} = draw_car(xx(t,1:2),[0,0,1],alphas(s));
            car_c1{s} = draw_car(xx(t,13:14),'g',alphas(s));
        end

        y = solution_set{1}{t}(14);
        axis([-30,+30,y-10,y+50]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(0.1);
        end
        for s = 1:L
%             delete(spot_a1{s});
            delete(car_a1{s});
%             delete(spot_b1{s});
%             delete(car_b1{s});
%             delete(spot_b2{s});
            delete(car_b2{s});
            delete(car_c1{s});
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
    xb2 = x(9:10);
    const = zeros(0,1);
%     const = [left_lane_constraint(xa1)+0.1];
end

function const = running_ineq_constraint_b1(x)
    xa1 = x(1:2);
    xb1 = x(5:6);
    xb2 = x(9:10);
    xc1 = x(13:14);
    const = [avoid_constraint(xa1,xb1)];
%              right_lane_ineq_constraint(xb1);
%              left_lane_ineq_constraint(xb1)];
end

function const = running_ineq_constraint_b2(x)
    xa1 = x(1:2);
    xb1 = x(5:6);
    xb2 = x(9:10);
    xc1 = x(13:14);
    const = [avoid_constraint(xa1,xb2)];
%     const = zeros(0,1);
%             right_lane_ineq_constraint(xb2);
%             left_lane_ineq_constraint(xb2)];
end

function const = running_ineq_constraint_c1(x)
    xa1 = x(1:2);
    xb1 = x(5:6);
    xb2 = x(9:10);
    xc1 = x(13:14);
    const = [avoid_constraint(xa1,xc1)];
end



%% Terminal ineq Constraints

function const = terminal_ineq_constraint_a1(x)
    const = [running_ineq_constraint_a1(x)];
end

function const = terminal_ineq_constraint_b1(x)    
    const = [running_ineq_constraint_b1(x)];
end

function const = terminal_ineq_constraint_b2(x)
    const = [running_ineq_constraint_b2(x)];
end

function const = terminal_ineq_constraint_c1(x)
    const = [running_ineq_constraint_c1(x)];
end

%% Terminal constraints 

function const = terminal_constraint_a1(x)
    xa1 = x(1:2);
    const = left_lane_constraint(xa1);
end

function const = terminal_constraint_b1(x)
    xa1 = x(1:2);
    xb1 = x(5:6);
    const = [left_lane_constraint(xb1)];
end

function const = terminal_constraint_b2(x)
    xb2 = x(9:10);
    const = left_lane_constraint(xb2);
end

function const = terminal_constraint_c1(x)
    xc1 = x(13:14);
    const = right_lane_constraint(xc1);
end

function val = terminal_cost(~)
    val = 0;
end

%% Running costs
function val = running_cost_a1(xu)
    global gv_a1 speed_pen lane_pen;
    global R;
    xa1 = xu(1:4);
    ua1 = xu(17:18);
    val = ua1'*R*ua1 + speed_pen*speed_cost(xa1,gv_a1);
end

function val = running_cost_a1_ll(xu)
    global gv_a1;
    global speed_pen;
    global R;
    xa1 = xu(1:4);
    ua1 = xu(17:18);
    val = ua1'*R*ua1 + speed_pen*speed_cost(xa1,gv_a1) + 100*left_lane_cost(xa1);
end

function val = running_cost_b1(xu)
    global gv_b1;
    global speed_pen;
    global R;
    xb1 = xu(5:8);
    ub1 = xu(19:20);
    val = ub1'*R*ub1 + speed_pen*speed_cost(xb1,gv_b1) + 100*left_lane_cost(xb1);
end

function val = running_cost_b2(xu)
    global gv_b2;
    global speed_pen;
    global R;
    xb2 = xu(9:12);
    ub2 = xu(21:22);
    val = ub2'*R*ub2 + speed_pen*speed_cost(xb2,gv_b2) + 100*left_lane_cost(xb2);
end

function val = running_cost_c1(xu)
    global gv_c1;
    global speed_pen;
    global R;
    xc1 = xu(13:16);
    uc1 = xu(23:24);
    val = uc1'*R*uc1 + speed_pen*speed_cost(xc1,gv_c1);
end


%% Polite costs
function val = running_cost_a1_polite(xu)
    global prob_b1;
    global prob_b2;
    pe = (1-prob_b1)*(1-prob_b2);
    pe1 = (1-prob_b2)*prob_b1;
    pe2 = (1-prob_b1)*prob_b2;
    p12 = prob_b1*prob_b2;
    
    r1 = prob_b1/(1-prob_b1);
    r2 = prob_b2/(1-prob_b2);
    
%     val = pe*running_cost_a1(xu) + ...
%         pe1*(running_cost_a1(xu)+running_cost_b1(xu)) + ...
%         pe2*(running_cost_a1(xu)+running_cost_b2(xu)) + ...
%         p12*(running_cost_b1(xu)+running_cost_b2(xu));
    
    
    val = running_cost_a1(xu) + r1*running_cost_b1(xu) + r2*running_cost_b2(xu);
end

function val = running_cost_a1_polite_ll(xu)
    global prob_b1;
    global prob_b2;
    pe = (1-prob_b1)*(1-prob_b2);
    pe1 = (1-prob_b2)*prob_b1;
    pe2 = (1-prob_b1)*prob_b2;
    p12 = prob_b1*prob_b2;
    
    r1 = prob_b1/(1-prob_b1);
    r2 = prob_b2/(1-prob_b2);
    val = running_cost_a1_ll(xu) + r1*running_cost_b1(xu) + r2*running_cost_b2(xu);
%     
%     val = pe*running_cost_a1_ll(xu) + ...
%         pe1*(running_cost_a1_ll(xu)+running_cost_b1(xu)) + ...
%         pe2*(running_cost_a1_ll(xu)+running_cost_b2(xu)) + ...
%         p12*(running_cost_b1(xu)+running_cost_b2(xu));
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
%     val = goal-x(4);
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

function const = avoid_constraint(x1,x2)
    global avoid_r;
    delt = x1(1:2)-x2(1:2);
    const = sqrt(delt'*delt)-avoid_r;
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

function pt = draw_car(center,color,alpha)
    global half_l half_w;
    pt = patch('XData',[center(1)-half_w,center(1)-half_w,center(1)+half_w,center(1)+half_w,center(1)-half_w],...
               'YData',[center(2)-half_l,center(2)+half_l,center(2)+half_l,center(2)-half_l,center(2)-half_l],...
               'FaceAlpha',alpha,'FaceColor', color);
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
%     plot([border_left/2, border_left/2],[-200,200],'--k');
    plot([0, 0],[-200,200],'--k');
%     plot([border_right/2, border_right/2],[-200,200],'--k');
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
    global dt;
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
    xc1 = xu(13:16);
    
    ua1 = xu(17:18);
    ub1 = xu(19:20);
    ub2 = xu(21:22);
    uc1 = xu(23:24);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    zb2 = xb2 + dt*unicycle_dyn(xb2,ub2);
    zc1 = xc1 + dt*unicycle_dyn(xc1,uc1);
    next = [za1;zb1;zb2;zc1];
end