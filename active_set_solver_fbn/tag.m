%% Tag
clear; clc; close all;

global dt;
dt = 0.25;

global linear_dyn;
linear_dyn = true;

global barrier_center barrier_radius;
barrier_center = [0,0];
barrier_radius = 3; 

global catch_weight;
catch_weight = .1;

global goal_a;
goal_a = [5;5];

global ped_r avoid_r;

ped_r = 0.5;
avoid_r = (ped_r+ped_r);


global T
T = 50;

x0 = [-5;1;0;0;
       5;1;0;0];
  

rollout_x0 = x0;
  
m{1} = 2;
m{2} = 2;

N = 2;
n = N*4;


for t = 1:T
    l{t,1} = @running_cost_a;
    l{t,2} = @running_cost_b;
  
    g{t,1} = @running_ineq_constraint_a;
    g{t,2} = @running_ineq_constraint_b;
    
%     g{t,1} = @empty_constraint;
%     g{t,2} = @empty_constraint;
    
    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;


end
l{T+1,1} = @empty_cost;
l{T+1,2} = @empty_cost;

g{T+1,1} = @terminal_ineq_constraint_a;
g{T+1,2} = @terminal_ineq_constraint_b;

% g{T+1,1} = @empty_constraint;
% g{T+1,2} = @empty_constraint;

% h{T+1,1} = @terminal_constraint_a;
% h{T+1,2} = @terminal_constraint_b;

h{T+1,1} = @empty_constraint;
h{T+1,2} = @empty_constraint;


%% Run 
params.use_ws = false;
params.use_initialization = false;
params.linear = false;
params.debug_plot = false;
params.open_loop = false;
[xval,ws,duration] = fb_game_solver_active_set(@f, h, g, l, n, m, N, T, x0, rollout_x0, params, @plot_solution);

plot_solution(xval,0);



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
    global T;
    draw_barrier();

    for t = 1:T+1
        xx(t,:) = xval{t};
        
        p1 = draw_ped(xx(t,1:2),'b');
        p2 = draw_ped(xx(t,5:6),'r');
      
        
        axis([-20,+20,-20,20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(.01);
        end
        delete(p1);
        delete(p2);
         
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
        writerObj = VideoWriter('drive_slow','MPEG-4');
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
        for s = 1:L
            xx(t,:) = solution_set{s}{t};    
            car_b1{s} = draw_car(xx(t,5:6),[1,0,0],alphas(s));
        end
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

        y = solution_set{1}{t}(6);
        axis([-30,+30,y-30,y+30]);
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
            delete(car_b1{s});
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

function const = running_ineq_constraint_a(x)
    xa = x(1:2);
    xb = x(5:6);
    const = [avoid_constraint(xa,xb); avoid_barrier(xa)];
end

function const = running_ineq_constraint_b(x)
    xa = x(1:2);
    xb = x(5:6);  
    const = zeros(0,1);
    const = avoid_barrier(xb);
end


%% Terminal ineq Constraints

function const = terminal_ineq_constraint_a(x)
    const = [running_ineq_constraint_a(x)];
end

function const = terminal_ineq_constraint_b(x)    
    const = [running_ineq_constraint_b(x)];
end

%% Terminal constraints 

function const = terminal_constraint_a(x)
    xa = x(1:2);
    global goal_a;
    const = xa-goal_a;
end

function const = terminal_constraint_b(x)
    xb = x(5:6);
    const = zeros(0,1);
end

function val = terminal_cost(~)
    val = 0;
end

%% Running costs
function val = running_cost_a(xu)
    global catch_weight;
    xa = xu(1:4);
    xb = xu(5:8);
    ua = xu(9:10);
    val = ua'*ua + 0.01*xa(1:2)'*xa(1:2);
end

function val = running_cost_b(xu)
    global catch_weight;
    xa = xu(1:4);
    xb = xu(5:8);
    ub = xu(11:12);
    val = ub'*ub + catch_weight*catch_cost(xa,xb);
end

%% Helpers

function const = avoid_constraint(x1,x2)
    global avoid_r;
    delt = x1(1:2)-x2(1:2);
    const = sqrt(delt'*delt)-avoid_r;
end

function const = avoid_barrier(x)
    global ped_r barrier_radius;
    delt = x(1:2);
    const = sqrt(delt'*delt)-(ped_r+barrier_radius);
end

function cost = catch_cost(x1,x2)
    delt = x1(1:2) - x2(1:2);
    cost = sqrt(delt'*delt);
end

function pt = draw_ped(center,color)
    global ped_r;
    pt = viscircles(center,ped_r,'Color', color, 'LineStyle', '-');
end

function draw_barrier()
    global barrier_center barrier_radius;
    barrier = viscircles(barrier_center,barrier_radius,'Color', 'k','LineStyle','-');
end

function vec = unicycle_dyn(x,u)
    global linear_dyn;
    vec = x(4)*cos(x(3));
    vec = [vec; x(3)*sin(x(3))];
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
    
    ua = xu(9:10);
    ub = xu(11:12);
    
    za = xa + dt*unicycle_dyn(xa,ua);
    zb = xb + dt*unicycle_dyn(xb,ub);
    
    next = [za;zb];
end