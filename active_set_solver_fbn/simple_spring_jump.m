%% Is there a car?
clear; clc; close all;

global dt;
dt = 0.1;

global linear_dyn;
linear_dyn = true;

global T
T = 40;

global l k1 k2 g m1 m2;

m1 = 1;
m2 = .5;

l = 4;
k1 = 50;
k2 = 1;
g = 9.81;

A = [0 0 1 0; 
    0 0  0 1;
    -k1/m1 k1/m1 -k2/m1 k2/m1;
    k1/m2 -k1/m2 k2/m2 -k2/m2];
c = [0; 0; l*k1/m1-g; -l*k1/m2-g];
Bu = [0; 0; 1/m1; -1/m2];
Bv = [0; 0; 0; 1/m2];
[V,J] = jordan(A);

A1 = J(1:2,1:2);
A2 = J(3:end,3:end);


Ad = expm(A*dt);
mid1 = [dt, -dt*dt/2; 0 dt];
mid2 = inv(-A2)*(expm(-dt*A2)-eye(2));
mid = blkdiag(mid1,mid2);

B_base = Ad*V*mid*inv(V);
Bud = B_base*Bu;
Bvd = B_base*Bv;
cd = B_base*c;

%% TEST DISCRETIZED SYSTEM
x0 = [5;0;0;0];

u = -100;
v = 9;
steps = T;
ddt = dt/steps;

x1 = Ad*x0 + cd+Bvd*v+Bud*u;

xx{1} = x0;
yy(:,1) = x0;
tot = 0;
for t = 1:steps
    tot = tot+ddt;
    xx{t+1} = xx{t}+ddt*(A*xx{t}+c+Bv*v+ Bu*u);
    yy(:,t+1) = xx{t+1};
end

plot_solution(xx);
%%
figure;
plot(1:(steps+1),yy(1,:));
hold on;
plot(1:(steps+1),yy(2,:));
plot(steps+1,real(x1(1)),'Marker','x','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
plot(steps+1,real(x1(2)),'Marker','x','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);

%%

At = eye(4) + dt*A;


rank([Bud Ad*Bud Ad*Ad*Bud Ad*Ad*Ad*Bud])
rank([Bu At*Bu At*At*Bu At*At*At*Bu])
%%



x0 = [7;3;0;0];  %vr
  

rollout_x0 = x0;

m{1} = 1;
m{2} = 1;
m{3} = 1;

N = 3;
n = 4;

full_scope{1} = true;
full_scope{2} = false;
full_scope{3} = false;

for t = 1:T
    
    lf{t,1} = @running_cost_robot;
    lf{t,2} = @running_cost_normal;
    lf{t,3} = @running_cost_limits;
  
    gf{t} = assign_running_ineq_constraints();
    gi{t} = assign_running_ineq_ownerships();
    
    h{t} = @empty_constraint;
    hi{t} = {};

end
lf{T+1,1} = @empty_cost;
lf{T+1,2} = @running_cost_normal;
lf{T+1,3} = @running_cost_limits;

gf{T+1} = assign_terminal_ineq_constraints();
gi{T+1} = assign_terminal_ineq_ownerships();

h{T+1} = @terminal_constraint;
hi{T+1}{1} = 1;

%% Run 

params.use_ws = false;
params.open_loop = false;
params.use_initialization = false;
params.debug_plot = true;
params.linear = false;
params.contact = true;
params.full_scope = full_scope;

[xval,ws,durs] = fb_game_solver_active_set_shared(@f, h, hi, gf, gi, lf, n, m, N, T, x0, rollout_controls, params, @plot_solution);
solution_set{1} = xval;
durations{1} = durs;

%% Plotting

function plot_solution(xval,x0)
    make_vid = false;
    close all;
    fh = draw_environment(); 
    hold on;
    if make_vid
        writerObj = VideoWriter('spring_jump','MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    
    global T m1 m2 dt;
    

    for t = 1:T+1
        xx(t,:) = xval{t};
        
        p1 = xx(t,1);
        p2 = xx(t,2);
        
        mass_1 = viscircles([t,p1],m1,'Color', 'b','LineStyle','-');
        mass_2 = viscircles([t,p2],m2,'Color', 'b','LineStyle','-');
        

        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        else
            pause(dt);
        end
        delete(mass_1);
        delete(mass_2);
        
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
    normal = xu(11:12);
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
    normal = xu(11:12);
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
    u = xu(9:10);
    normal = xu(11:12);
    tangential = xu(13:14);
    const = [0.001 + mu*normal(1) - tangential(1);
             0.001 + mu*normal(1) + tangential(1)];
%              0.001 - mu*normal(1) - tangential(1);
%              0.001 - mu*normal(1) + tangential(1)];
end

function const = running_ineq_constraint_tangential_right(xu)
    global mu;
    x = xu(1:8);
    u = xu(9:10);
    normal = xu(11:12);
    tangential = xu(13:14);
    const = [0.001 + mu*normal(2) - tangential(2);
             0.001 + mu*normal(2) + tangential(2)];
%              0.001 - mu*normal(2) - tangential(2);
%              0.001 - mu*normal(2) + tangential(2)];
end


function region = running_ineq_region_tangential_left(xu)
    normal_left = xu(11);
    region = [1;1];
%     region = [normal_left; normal_left; -normal_left-1e-6; -normal_left-1e-6];
end

function region = running_ineq_region_tangential_right(xu)
    normal_right = xu(12);
    region = [1;1];
%     region = [normal_right; normal_right; -normal_right-1e-6; -normal_right-1e-6];
end

function const = running_ineq_constraint_limits_left(xu)
    x = xu(1:8);
    limits = xu(17:18);
    const = limits(1);
end

function region = running_ineq_region_limits_left(xu)
    region = 1;
end

function const = running_ineq_constraint_limits_right(xu)
    x = xu(1:8);
    limits = xu(17:18);
    const = limits(2);
end

function region = running_ineq_region_limits_right(xu)
    region = 1;
end


function const = running_ineq_constraint(xu)
    const = [running_ineq_constraint_robot(xu);
            running_ineq_constraint_normal_left(xu);
            running_ineq_constraint_normal_right(xu);
            running_ineq_constraint_tangential_left(xu);
            running_ineq_constraint_tangential_right(xu);
            running_ineq_constraint_limits_left(xu);
            running_ineq_constraint_limits_right(xu)];
end

function region = running_ineq_region(xu)
    region = [running_ineq_region_robot(xu);
        running_ineq_region_normal_left(xu);
        running_ineq_region_normal_right(xu);
        running_ineq_region_tangential_left(xu);
        running_ineq_region_tangential_right(xu);
        running_ineq_region_limits_left(xu);
        running_ineq_region_limits_right(xu)];     
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
    region = zeros(0,1);
end

function region = terminal_ineq_region_normal_right(x)
    region = zeros(0,1);
end

function region = terminal_ineq_region_tangential_left(x)
    region = zeros(0,1);
end

function region = terminal_ineq_region_tangential_right(x)
    region = zeros(0,1);
end

function const = terminal_ineq_constraint(x)
    const= [terminal_ineq_constraint_robot(x);
            terminal_ineq_constraint_normal_left(x);
            terminal_ineq_constraint_normal_right(x);
            terminal_ineq_constraint_tangential_right(x);
            terminal_ineq_constraint_tangential_left(x)];
end

function region = terminal_ineq_region(x)
    region = [terminal_ineq_region_robot(x);
            terminal_ineq_region_normal_left(x);
            terminal_ineq_region_normal_right(x);
            terminal_ineq_region_tangential_right(x);
            terminal_ineq_region_tangential_left(x)];
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

function const = terminal_constraint(x)
    const = [terminal_constraint_robot(x);
            terminal_constraint_normal_left(x);
            terminal_constraint_normal_right(x);
            terminal_constraint_tangential_left(x);
            terminal_constraint_tangential_right(x)];
end

%% Running costs
function val = running_cost_robot(xu)
    x = xu(1:8);
    u = xu(9:10);
    val = 0.5*u'*u;
end

function val = running_cost_cheat(xu)
    x = xu(1:8);
    v = xu(15:16);
    val = 0.5*v'*v;
end

function val = running_cost_limits_left(xu)
    global l;
    x = xu(1:8);
    val = 0.5*(x(5)+l)^2;
end

function val = running_cost_limits_right(xu)
    global l;
    x = xu(1:8);
    val = 0.5*(x(7)+l)^2;
end

function val = running_cost_normal_left(xu)
    global wall_angle;
    x = xu(1:8);
%     u = xu(9:10);
%     normal_left = xu(11);
%     tangential = xu(13:14);
    point_left = [x(1); x(2)] + x(5) *  [-sin(wall_angle); -cos(wall_angle)];
    val = 0.5*(point_left'*[cos(wall_angle);  sin(wall_angle)])^2;
%     val = 0.5*normal_left'*normal_left;
end

function val = running_cost_normal_right(xu)
    global wall_angle;
    x = xu(1:8);
%     u = xu(9:10);
%     normal_right = xu(12);
%     tangential = xu(13:14);
    point_right = [x(1); x(2)] + x(7) * [ sin(wall_angle); -cos(wall_angle)];
    val = 0.5*(point_right'*[cos(wall_angle);  sin(-wall_angle)])^2;
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

function val = terminal_cost_limits_left(xu)
    global l;
    x = xu(1:8);
    val = 0.5*(x(5)+l)^2;
end

function val = terminal_cost_limits_right(xu)
    global l;
    x = xu(1:8);
    val = 0.5*(x(7)+l)^2;
end


%% Helpers

function fh = draw_environment()
    fh = figure;
    global T
    hold on;
    plot([-1, T+1],[0,0],'-k','linewidth', 4);
    axis('equal');
    axis([-1,T+1,-1,T+1]);
end


function next = f(xu)
    global dt wall_angle Mc Me l C K g;
    X = xu(1:8);
    U = xu(9:10);
    Y = xu(11:12);
    Z = xu(13:14);
    V = xu(15:16);
    lims = xu(17:18);
                                
    left_spring_force = K*(X(5)-l) + C*X(6)+U(1);
    right_spring_force = K*(X(7)-l) + C*X(8)+U(2);
    
    left_com_accel = -[sin(wall_angle); cos(wall_angle)] * left_spring_force/(Mc+Me) + [cos(wall_angle); -sin(wall_angle)] * Z(1) / (Mc+2*Me);
    right_com_accel = -[-sin(wall_angle); cos(wall_angle)] * right_spring_force/(Mc+Me) + [cos(wall_angle); sin(wall_angle)] * Z(2) / (Mc+2*Me);
    
    left_length_accel = (-Y(1) - left_spring_force + lims(1)) / Me;
    right_length_accel = (-Y(2) - right_spring_force + lims(2)) / Me;
    
    com_accel = left_com_accel+right_com_accel+[V(1);V(2)-g];
    
    next = X + dt * [X(3)+dt/2*com_accel(1);
                     X(4)+dt/2*com_accel(2);
                     com_accel;
                     X(6)+dt/2*left_length_accel;
                     left_length_accel;
                     X(8)+dt/2*right_length_accel;
                     right_length_accel];
end