%% Is there a car?
clear; close all;

global dt;
dt = 0.04;

global linear_dyn;
linear_dyn = true;

global T
T = 50;

global wall_angle Mc Me bar mu l K C g h0;

wall_angle = pi/4;

Mc = 1;
Me = 0.1;
bar = 6;

mu = 0.1;
l = 4;
K = 100;
C = 1;

g = 9.81;
h0 = 15;

x0 = [7.5; %cx
    h0;   %cy
    0;   %vx
    -3;   %vy
    l;   %dl
    0;   %vl
    l;   %dr
    0];  %vr
  

rollout_x0 = x0;

m{1} = 1;
m{2} = 1;
m{3} = 1;
m{4} = 1;
m{5} = 1;
N = 5;
n = 8;

for t = 1:T
    
    lf{t,1} = @running_cost_robot;
    lf{t,2} = @running_cost_normal_left;
    lf{t,3} = @running_cost_normal_right;
    lf{t,4} = @running_cost_tangential_left;
    lf{t,5} = @running_cost_tangential_right;
  
    gf{t,1} = @running_ineq_constraint_robot;
    gf{t,2} = @running_ineq_constraint_normal_left;
    gf{t,3} = @running_ineq_constraint_normal_right;
    gf{t,4} = @running_ineq_constraint_tangential_left;
    gf{t,5} = @running_ineq_constraint_tangential_right;
    
    gr{t,1} = @running_ineq_region_robot;
    gr{t,2} = @running_ineq_region_normal_left;
    gr{t,3} = @running_ineq_region_normal_right;
    gr{t,4} = @running_ineq_region_tangential_left;
    gr{t,5} = @running_ineq_region_tangential_right;
    
    h{t,1} = @empty_constraint;
    h{t,2} = @empty_constraint;
    h{t,3} = @empty_constraint;
    h{t,4} = @empty_constraint;
    h{t,5} = @empty_constraint;
    
end
lf{T+1,1} = @empty_cost;
lf{T+1,2} = @running_cost_normal_left;
lf{T+1,3} = @running_cost_normal_right;
lf{T+1,4} = @terminal_cost_tangential_left;
lf{T+1,5} = @terminal_cost_tangential_right;

gf{T+1,1} = @terminal_ineq_constraint_robot;
gf{T+1,2} = @terminal_ineq_constraint_normal_left;
gf{T+1,3} = @terminal_ineq_constraint_normal_right;
gf{T+1,4} = @terminal_ineq_constraint_tangential_left;
gf{T+1,5} = @terminal_ineq_constraint_tangential_right;

gr{T+1,1} = @terminal_ineq_region_robot;
gr{T+1,2} = @terminal_ineq_region_normal_left;
gr{T+1,3} = @terminal_ineq_region_normal_right;
gr{T+1,4} = @terminal_ineq_region_tangential_left;
gr{T+1,5} = @terminal_ineq_region_tangential_right;

h{T+1,1} = @terminal_constraint_robot;
h{T+1,2} = @terminal_constraint_normal_left;
h{T+1,3} = @terminal_constraint_normal_right;
h{T+1,4} = @terminal_constraint_tangential_left;
h{T+1,5} = @terminal_constraint_tangential_right;


%% Run 

params.use_ws = false;
params.open_loop = false;
params.use_initialization = false;
params.debug_plot = true;
params.linear = false;
params.contact = true;


[xval,ws,durs] = fb_game_solver_active_set_mod(@f, h, gf, gr, lf, n, m, N, T, x0, rollout_x0, params, @plot_solution);
solution_set{1} = xval;
durations{1} = durs;

%% Plotting

function plot_solution(xval,x0)
    make_vid = false;
    close all;
    fh = draw_environment(); 
    hold on;
    if make_vid
        writerObj = VideoWriter('wallball','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    
    global T wall_angle;
    

    for t = 1:T+1
        xx(t,:) = xval{t};
        
        px = xx(t,1);
        py = xx(t,2);
        dl = xx(t,5);
        dr = xx(t,7);
        
        pp = fimplicit(@(zx,zy) ((zx-px)*cos(-wall_angle) + (zy-py)*sin(-wall_angle))^2 / (dr^2) + ((zx-px)*sin(-wall_angle) - (zy-py)*cos(-wall_angle))^2 / (dl^2) - 1, '-b');

        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(.01);
        end
        delete(pp);
        
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
        writerObj = VideoWriter('merge_in_front','MPEG-4');
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

function const = running_ineq_constraint_robot(x)
    global l bar;
    x = x(1:8);
%     u = x(9:10);
%     normal = x(11:12);
%     tangential = x(13:14);
    const = zeros(0,1);
    const = x(2) - l - bar;
end

function region = running_ineq_region_robot(x)
    region = 1;
end

function const = running_ineq_constraint_normal_left(xu)
    x = xu(1:8);
    normal = xu(10:11);
%     global wall_angle;
%     
%     point_left = [x(1); x(2)] + x(5) *  [-sin(wall_angle); -cos(wall_angle)];
%     const = [point_left(2)  - tan(-wall_angle)*point_left(1)];
    const = normal(1);
end

function region = running_ineq_region_normal_left(xu)
    region = 1;
end

function const = running_ineq_constraint_normal_right(xu)
    x = xu(1:8);
    normal = xu(10:11);
%     global wall_angle;
%     
%     point_right = [x(1); x(2)] + x(7) * [ sin(wall_angle); -cos(wall_angle)];
%     const = [point_right(2) - tan(wall_angle)*point_right(1)];
    const = normal(2);
end

function region = running_ineq_region_normal_right(xu)
    region = 1;
end

function const = running_ineq_constraint_tangential_left(xu)
    global mu;
    x = xu(1:8);
    u = xu(9);
    normal = xu(10:11);
    tangential = xu(12:13);
    const = [0.001 + mu*normal(1) - tangential(1);
             0.001 + mu*normal(1) + tangential(1)];
%              0.001 - mu*normal(1) - tangential(1);
%              0.001 - mu*normal(1) + tangential(1)];
end

function const = running_ineq_constraint_tangential_right(xu)
    global mu;
    x = xu(1:8);
    u = xu(9);
    normal = xu(10:11);
    tangential = xu(12:13);
    const = [0.001 + mu*normal(2) - tangential(2);
             0.001 + mu*normal(2) + tangential(2)];
%              0.001 - mu*normal(2) - tangential(2);
%              0.001 - mu*normal(2) + tangential(2)];
end


function region = running_ineq_region_tangential_left(xu)
    normal_left = xu(10);
    region = [1;1];
%     region = [normal_left; normal_left; -normal_left-1e-6; -normal_left-1e-6];
end

function region = running_ineq_region_tangential_right(xu)
    normal_right = xu(11);
    region = [1;1];
%     region = [normal_right; normal_right; -normal_right-1e-6; -normal_right-1e-6];
end



%% Terminal ineq Constraints

function const = terminal_ineq_constraint_robot(x)
    const = [running_ineq_constraint_robot(x)];
end

function const = terminal_ineq_constraint_normal_left(x)    
%     const = [running_ineq_constraint_normal_left(x)];
    const = zeros(0,1);
end

function const = terminal_ineq_constraint_normal_right(x)    
%     const = [running_ineq_constraint_normal_right(x)];
    const = zeros(0,1);
end

function const = terminal_ineq_constraint_tangential_left(~)
    const = zeros(0,1);
end
function const = terminal_ineq_constraint_tangential_right(~)
    const = zeros(0,1);
end

function region = terminal_ineq_region_robot(x)
    region = 1;
end

function region = terminal_ineq_region_normal_left(x)
    region = 1;
end

function region = terminal_ineq_region_normal_right(x)
    region = 1;
end

function region = terminal_ineq_region_tangential_left(x)
    region = zeros(0,1);
end

function region = terminal_ineq_region_tangential_right(x)
    region = zeros(0,1);
end

%% Terminal constraints 

function const = terminal_constraint_robot(x)
    const = zeros(0,1);
end

function const = terminal_constraint_normal_left(x)
    const = zeros(0,1);
end

function const = terminal_constraint_normal_right(x)
    const = zeros(0,1);
end


function const = terminal_constraint_tangential_left(x)
    const = zeros(0,1);
end

function const = terminal_constraint_tangential_right(x)
    const = zeros(0,1);
end

function val = terminal_cost(~)
    val = 0;
end

%% Running costs
function val = running_cost_robot(xu)
    x = xu(1:8);
    u = xu(9);
    val = 0.5*u*u;
end

function val = running_cost_normal_left(xu)
    global wall_angle;
    x = xu(1:8);
%     u = xu(9:10);
%     normal_left = xu(11);
%     tangential = xu(13:14);
    point_left = [x(1); x(2)] + x(5) *  [-sin(wall_angle); -cos(wall_angle)];
    val = 10*0.5*(point_left'*[cos(wall_angle);  sin(wall_angle)])^2;
%     val = 0.5*normal_left'*normal_left;
end

function val = running_cost_normal_right(xu)
    global wall_angle;
    x = xu(1:8);
%     u = xu(9:10);
%     normal_right = xu(12);
%     tangential = xu(13:14);
    point_right = [x(1); x(2)] + x(7) * [ sin(wall_angle); -cos(wall_angle)];
    val = 10*0.5*(point_right'*[cos(wall_angle);  sin(-wall_angle)])^2;
%     val = 0.5*normal_right'*normal_right;
end

function val = running_cost_tangential_left(xu)
    global wall_angle;
    x = xu(1:8);
%     u = xu(9:10);
%     normal = xu(11:12);
%     tangential = xu(13:14);
    val = 0.5*([x(3);x(4)]'*[cos(wall_angle);  sin(-wall_angle)])^2;
end

function val = running_cost_tangential_right(xu)
    global wall_angle;
    x = xu(1:8);
%     u = xu(9:10);
%     normal = xu(11:12);
%     tangential = xu(13:14);
    val = 0.5*([x(3);x(4)]'*[cos(wall_angle);  sin(wall_angle)])^2;
end

function val = terminal_cost_tangential_left(xu)
    global wall_angle;
    x = xu(1:8);
    val = 0.5*([x(3);x(4)]'*[cos(wall_angle);  sin(-wall_angle)])^2;
end

function val = terminal_cost_tangential_right(xu)
    global wall_angle;
    x = xu(1:8);
    val = 0.5*([x(3);x(4)]'*[cos(wall_angle);  sin(wall_angle)])^2;
end


%% Helpers

function fh = draw_environment()
    global bar wall_angle l;
    fh = figure;
    hold on;
    bound = 15;
    bound = bound + l + 3;
    plot([-bound,bound],[bar, bar], '--r');
    h=fill([-bound -bound bound bound],[-bound bar bar -bound],'red');
    h.FaceAlpha=0.5;
    plot([-bound,0],tan(-wall_angle)*[-bound,0],'-k','linewidth', 4)
    plot([0,bound],tan(wall_angle)*[0,bound],'-k','linewidth', 4)
    axis('equal')
    axis([-bound,bound,-bound/2,3/2*bound]);
end


function next = f(xu)
    global dt wall_angle Mc Me l C K g;
    X = xu(1:8);
    U = xu(9);
    Y = xu(10:11);
    Z = xu(12:13);
                                
    left_spring_force = K*(X(5)-l) + C*X(6);
    right_spring_force = K*(X(7)-l) + C*X(8);
    
    left_com_accel = -[sin(wall_angle); cos(wall_angle)] * left_spring_force/(Mc+Me) + [cos(wall_angle); -sin(wall_angle)] * Z(1) / (Mc+2*Me);
    right_com_accel = -[-sin(wall_angle); cos(wall_angle)] * right_spring_force/(Mc+Me) + [cos(wall_angle); sin(wall_angle)] * Z(2) / (Mc+2*Me);
    
    left_length_accel = (-Y(1) - left_spring_force) / Me;
    right_length_accel = (-Y(2) - right_spring_force) / Me;
    
    com_accel = left_com_accel+right_com_accel+[0;-g+U];
    
    next = X + dt * [X(3)+dt/2*com_accel(1);
                     X(4)+dt/2*com_accel(2);
                     com_accel;
                     X(6)+dt*left_length_accel;
                     left_length_accel;
                     X(8)+dt*right_length_accel;
                     right_length_accel];
end