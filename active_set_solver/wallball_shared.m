%% Is there a car?
clear; clc; close all;

global dt;
dt = 0.1;

global linear_dyn;
linear_dyn = true;

global T
T = 50;

global spring rc_on;
spring = true;
rc_on = false;

global wall_angle Mc Me bar mu l K C g h0;

wall_angle = pi/4;

Mc = 5;
Me = 0.1;
bar = 6;

mu = 0.3;
l = 4;
K = 100;
C = 1;

g = 9.81;
h0 = 15;

x0 = [7.5; %cx
    h0;   %cy
    2;   %vx
    -3;   %vy
    l;   %dl
    0;   %vl
    l;   %dr
    0];  %vr
  

if spring
    rollout_controls{1} = [0; 0];
else
    rollout_controls{1} = 15;
end
rollout_controls{2} = 0;
rollout_controls{3} = 0;
rollout_controls{4} = 0;
rollout_controls{5} = 0;

full_scope{1} = true;
full_scope{2} = false;
full_scope{3} = false;
full_scope{4} = false;
full_scope{5} = false;


if spring
    m{1} = 2;
else
    m{1} = 1;
end
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
      
    gf{t} = assign_running_ineq_constraint();
    gi{t} = assign_running_ineq_ownerships();
    
    h{t} = @empty_constraint;
    hi{t} = {};
    
end
lf{T+1,1} = @terminal_robot_cost;
lf{T+1,2} = @running_cost_normal_left;
lf{T+1,3} = @running_cost_normal_right;
lf{T+1,4} = @running_cost_tangential_left;
lf{T+1,5} = @running_cost_tangential_right;


gf{T+1} = assign_terminal_ineq_constraint();
gi{T+1} = assign_terminal_ineq_ownerships();

if ~rc_on
    gf{T+1} = {};
    gi{T+1} = {};
end


h{T+1} = @empty_constraint;
hi{T+1} = {};


%% Run 

params.use_ws = false;
params.open_loop = false;
params.use_initialization = false;
params.debug_plot = true;
params.linear = false;
params.full_scope = full_scope;
params.print_frequency = 50;
params.rollout_controls = rollout_controls;
params.wrap_val = inf;
mpc_iters = 1;

[xvals,durs,iters] = mpgp_fb_game_solver_active_set_shared(@f, h, hi, gf, gi, lf, n, m, N, T, mpc_iters, x0, params, @plot_solution);


solution_set{1} = xvals;
durations{1} = durs;
%%
plot_solution(xvals{1},0)

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
    const = zeros(0,1);
    const = x(2) - l - bar;
end

function const = running_ineq_constraint_normal_left(xu)
    x = xu(1:8);
    global spring
    if spring
        normal = xu(11:12);
    else
        normal = xu(10:11);
    end
    const = normal(1);
end

function const = running_ineq_constraint_normal_right(xu)
    x = xu(1:8);
    global spring
    if spring
        normal = xu(11:12);
    else
        normal = xu(10:11);
    end
    const = normal(2);
end

function const = running_ineq_constraint_tangential_left_1(xu)
    global mu;
    global spring;
    if spring
        x = xu(1:8);
        u = xu(9:10);
        normal = xu(11:12);
        tangential = xu(13:14);
    else
        x = xu(1:8);
        u = xu(9);
        normal = xu(10:11);
        tangential = xu(12:13);
    end
    const = [0.001 + mu*normal(1) - tangential(1)];
end

function const = running_ineq_constraint_tangential_left_2(xu)
    global mu;
    global spring;
    if spring
        x = xu(1:8);
        u = xu(9:10);
        normal = xu(11:12);
        tangential = xu(13:14);
    else
        x = xu(1:8);
        u = xu(9);
        normal = xu(10:11);
        tangential = xu(12:13);
    end
    const = [0.001 + mu*normal(1) + tangential(1)];
end

function const = running_ineq_constraint_tangential_right_1(xu)
    global mu;
    global spring;
    if spring
        x = xu(1:8);
        u = xu(9:10);
        normal = xu(11:12);
        tangential = xu(13:14);
    else
        x = xu(1:8);
        u = xu(9);
        normal = xu(10:11);
        tangential = xu(12:13);
    end
    const = [0.001 + mu*normal(2) - tangential(2)];
end

function const = running_ineq_constraint_tangential_right_2(xu)
    global mu;
    global spring;
    if spring
        x = xu(1:8);
        u = xu(9:10);
        normal = xu(11:12);
        tangential = xu(13:14);
    else
        x = xu(1:8);
        u = xu(9);
        normal = xu(10:11);
        tangential = xu(12:13);
    end
    const = [0.001 + mu*normal(2) + tangential(2)];
end


function const = assign_running_ineq_constraint()
    global rc_on;
    const{1} = @running_ineq_constraint_normal_left;
    const{2} = @running_ineq_constraint_normal_right;
    const{3} = @running_ineq_constraint_tangential_left_1;
    const{4} = @running_ineq_constraint_tangential_left_2;
    const{5} = @running_ineq_constraint_tangential_right_1;
    const{6} = @running_ineq_constraint_tangential_right_2;
    if rc_on
        const{7} = @running_ineq_constraint_robot;
    end
end

function ownerships = assign_running_ineq_ownerships()
    global rc_on;
    ownerships{1} = 2;
    ownerships{2} = 3;
    ownerships{3} = 4;
    ownerships{4} = 4;
    ownerships{5} = 5;
    ownerships{6} = 5;
    if rc_on
        ownerships{7} = 1;
    end
end

%% Terminal ineq Constraints

function const = terminal_ineq_constraint_robot(x)
    const = running_ineq_constraint_robot(x);
end

function const = assign_terminal_ineq_constraint()
    const{1} = @terminal_ineq_constraint_robot; 
end

function ownership = assign_terminal_ineq_ownerships(x)
    ownership{1} = 1;
end


%% Running costs
function val = running_cost_robot(xu)
    global spring;
    if spring
        x = xu(1:8);
        u = xu(9:10);
        normal = xu(11:12);
        tangential = xu(13:14);
    else
        x = xu(1:8);
        u = xu(9);
        normal = xu(10:11);
        tangential = xu(12:13);
    end
    val = 0.5*u'*u;
end

function val = terminal_robot_cost(x)
    global bar
    val = 100*(x(2)-2*bar)^2;
end

function val = running_cost_normal_left(xu)
    global wall_angle;
    x = xu(1:8);
    point_left = [x(1); x(2)] + x(5) *  [-sin(wall_angle); -cos(wall_angle)];
    surface_dist = point_left'*[cos(wall_angle);  sin(wall_angle)];
    val = 0.5*surface_dist^2;
end

function val = running_cost_normal_right(xu)
    global wall_angle;
    x = xu(1:8);
    point_right = [x(1); x(2)] + x(7) * [ sin(wall_angle); -cos(wall_angle)];
    surface_dist = point_right'*[cos(wall_angle);  sin(-wall_angle)];
    val = 0.5*surface_dist^2;
end

function val = running_cost_tangential_left(xu)
    global wall_angle;
    x = xu(1:8);
    tangential_vel = [x(3);x(4)]'*[cos(wall_angle);  sin(-wall_angle)];
    val = 0.5*tangential_vel^2;
end

function val = running_cost_tangential_right(xu)
    global wall_angle;
    x = xu(1:8);
    tangential_vel = [x(3);x(4)]'*[cos(wall_angle);  sin(wall_angle)];
    val = 0.5*tangential_vel^2;
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
    global dt wall_angle Mc Me l C K g spring;
    if spring
        X = xu(1:8);
        U = xu(9:10);
        Y = xu(11:12);
        Z = xu(13:14);
    else
        X = xu(1:8);
        U = xu(9);
        Y = xu(10:11);
        Z = xu(12:13);
    end
                                
%     left_spring_force = K*(X(5)-l) + C*X(6) + U1;
%     right_spring_force = K*(X(7)-l) + C*X(8) + U2;
%     
%     left_com_accel = -[sin(wall_angle); cos(wall_angle)] * left_spring_force/(Mc+Me) + [cos(wall_angle); -sin(wall_angle)] * Z(1) / (Mc+2*Me);
%     right_com_accel = -[-sin(wall_angle); cos(wall_angle)] * right_spring_force/(Mc+Me) + [cos(wall_angle); sin(wall_angle)] * Z(2) / (Mc+2*Me);
%     
%     left_length_accel = (-Y(1) - left_spring_force) / Me;
%     right_length_accel = (-Y(2) - right_spring_force) / Me;
%     
%     com_accel = left_com_accel+right_com_accel+[0;-g+U];
%     % forward euler
%     next = X + dt * [X(3)+dt/2*com_accel(1);
%                      X(4)+dt/2*com_accel(2);
%                      com_accel;
%                      X(6)+dt/2*left_length_accel;
%                      left_length_accel;
%                      X(8)+dt/2*right_length_accel;
%                      right_length_accel];
                 
    d35 = -K*sin(wall_angle)/(Mc+Me);
    d36 = -C*sin(wall_angle)/(Mc+Me);
    d37 = K*sin(wall_angle)/(Mc+Me);
    d38 = C*sin(wall_angle)/(Mc+Me);
    
    d45 = -K*cos(wall_angle)/(Mc+Me);
    d46 = -C*cos(wall_angle)/(Mc+Me);
    d47 = -K*cos(wall_angle)/(Mc+Me);
    d48 = -C*cos(wall_angle)/(Mc+Me);
    
    A = [0 0 1 0 0 0 0 0;
         0 0 0 1 0 0 0 0;
         0 0 0 0 d35 d36 d37 d38;
         0 0 0 0 d45 d46 d47 d48;
         0 0 0 0 0 1 0 0;
         0 0 0 0 -K/Me -C/Me 0 0;
         0 0 0 0 0 0 0 1;
         0 0 0 0 0 0 -K/Me -C/Me];
    c = [0; 0; 0; 2*K*l*cos(wall_angle)/(Me+Mc)-g; 0; K*l/Me; 0; K*l/Me];
    if spring
        BU = [0 0;
              0 0; 
              -sin(wall_angle)/(Me+Mc)  sin(wall_angle)/(Me+Mc); 
              -cos(wall_angle)/(Me+Mc)  -cos(wall_angle)/(Me+Mc); 
              0 0;
              -1/Me 0;
              0 0;
              0 -1/Me];
    else
        BU = [0;0;0;1;0;0;0;0];
    end
    BY = [0 0;
          0 0;
          0 0;
          0 0;
          0 0;
          -1/Me 0;
          0 0;
          0 -1/Me];
    BZ = [0 0;
          0 0;
          cos(wall_angle) cos(wall_angle);
          -sin(wall_angle) sin(wall_angle);
          0 0;
          0 0;
          0 0;
          0 0] / (Mc+2*Me);
      
    allB = [BU BY BZ c];
    sys = ss(A,allB, eye(8),zeros(8,size(allB,2)));
    
    sysd = c2d(sys,dt);
    next = sysd.A*X + sysd.B*[U;Y;Z;1];
    
%     next = expm(A*dt)*(X+BU*U + BY*Y + BZ*Z + c);
%     next = (eye(8)-dt*A)\(X + dt*(BU*U + BY*Y + BZ*Z + c));
%     disp('check');
    % backward euler
%     next  = (I-dt*A)\ (x + dt*(Bu + c))
end