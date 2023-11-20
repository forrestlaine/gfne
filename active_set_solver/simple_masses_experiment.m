%% Is there a car?
clear; clc; close all;

global dt;
dt = 1;

global linear_dyn;
linear_dyn = true;


global R;
global alpha;


R = 1;
alpha = 1;

global gv gw;
gv = 1/4;
gw = 1/2;

global T
T = 13;
mpc_iters = 1;

x0 = [5/4; gv; 0; gw];


rollout_controls{1} = [10];
rollout_controls{2} = [0];
m{1} = 1;
m{2} = 1;

full_scope{1} = true;
full_scope{2} = true;

N = 2;
n = N*2;

for t = 1:T
    lf{t,1} = @running_cost_a;
    lf{t,2} = @running_cost_b;
  
    gf{t} = assign_running_ineq_constraints();
    gi{t} = assign_running_ineq_ownerships();
    
    h{t} = @empty_constraint;
    hi{t} = {};
end

% lf{T+1,1} = @terminal_cost_a;
% lf{T+1,2} = @terminal_cost_b;
% 
% gf{T+1} = assign_running_ineq_constraints();
% gi{T+1} = assign_running_ineq_ownerships();
% 
% h{T+1} = @empty_constraint;
% hi{T+1} = {};


%% Run 

params.use_ws = false;
params.open_loop = false;
params.use_initialization = false;
params.debug_plot = true;
params.linear = false;
params.full_scope = full_scope;
params.print_frequency = 10;
params.rollout_controls = rollout_controls;
params.wrap_val = 36;

[xvals,durs,iters] = mpgp_fb_game_solver_active_set_shared(@f, h, hi, gf, gi, lf, n, m, N, T, mpc_iters, x0, params, @plot_solution);
solution_set{1} = xvals;
durations{1} = durs;

%%
S = numel(xvals);
for k = 1:S
    xval{k} = xvals{k}{1};
end

for k = S+2:S+T+1
    xval{k-1} = xvals{S}{k-S};
end

plot_solution(xval,0)
sol = 0;
for t = 1:numel(xval)
    sol = sol + norm(xval{t});
%     gf{t}{1}([xval{t};0;0])
end
sol
% draw_borders = true;
% plot_solution(xval,0);
%%

% load('can_merge');
% plot_solution(xval,0);


%% Plotting

function plot_solution(xval,x0)

    close all;
    fh = figure;
    hold on;
    
    T = numel(xval);
    for t = 1:T
        xx(t,:) = xval{t};
    end
        
    color_a = 'b';
    color_b = 'g';

    plot([1:T],xx(:,1),'-b');
    plot([1:T],xx(:,3),'-g');
end

%%

function const = empty_constraint(~)
    const = zeros(0,1);
end

function val = empty_cost(~)
    val = 0;
end

%% Running Constraints


function const = running_ineq_constraint_b(x)
    px = x(1);
    py = x(3);
    const = px-(py+1);
end

function const = assign_running_ineq_constraints()
    const{1} = @running_ineq_constraint_b;
end

function ownerships = assign_running_ineq_ownerships()
    ownerships{1} = [2]; 
end

%% Terminal ineq Constraints

function const = assign_terminal_ineq_constraints()
    const = assign_running_ineq_constraints();
end

function ownership = assign_terminal_ineq_ownerships()
    ownership = assign_running_ineq_ownerships();
end

%% Terminal constraints 

function val = terminal_cost(~)
    val = 0;
end

function val = terminal_cost_a(xu)
    global alpha gv;
    xa = xu(1:2);
    val = 0.5*(xa(2)-gv)^2 + alpha*terminal_cost_b(xu);
end

function val = terminal_cost_b(xu)
    global gw;
    xb = xu(3:4);
    val = 0.5*(xb(2)-gw)^2;
end

%% Running costs
function val = running_cost_a(xu)
    global R;
    global gv;
    global alpha;
    xa = xu(1:2);
    ua = xu(5);
    val = 0.5*ua'*R*ua + 0.5*(xa(2)-gv)^2 + alpha*running_cost_b(xu);
end

function val = running_cost_b(xu)
    global R;
    global gw;
    xb = xu(3:4);
    ub = xu(6);
    val = 0.5*ub'*R*ub + 0.5*(xb(2)-gw)^2;
end



%% Helpers


function vec = dyn(x,u)
    global dt;
    vec = [x(2)+dt/2*u; u];

end

function next = f(xu)
    global dt;
    xa = xu(1:2);
    xb = xu(3:4);
    
    ua = xu(5);
    ub = xu(6);
    
    za = xa + dt*dyn(xa,ua);
    zb = xb + dt*dyn(xb,ub);
    next = [za;zb];
end