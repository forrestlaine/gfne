clear; clc; close all;
import casadi.*

% UNIFORM ROD FALLING!

inf_const = 999999999999.9;

l = 1;
M = 1;
g = 9.81;
I = 1/12 * M*l*l;
mu = 0;

dt = 0.025;
T = 20;
num_iterations = 1;

X{1} = MX.sym(['X_' num2str(1)], 3);


%     y   dy   theta dtheta
x0 = [1.0; 0; -pi/4];

f_x = 0;

h_x = [X{1}-x0];
g_x = [];
vars_x = X{1};

states = X{1};
controls = [];
lams_l = [];
gams_l = [];
lams_r = [];
gams_r = [];

for t = 1:T
    % Introduce symbolic vars.
    U{t} = MX.sym(['U_' num2str(t)], 1);
    % normal forces [l; r] (positive is up)
    YL{t} = MX.sym(['YL_' num2str(t)], 1);
    YR{t} = MX.sym(['YR_' num2str(t)], 1);
    
    X{t+1} = MX.sym(['X_' num2str(t+1)], 3);

    vars_x = [vars_x; U{t}; X{t+1}];
    vars_y_l{t} = [YL{t}; X{t+1}];
    vars_y_r{t} = [YR{t}; X{t+1}];
   
    % Place vars in lists for easy access later.
    states = [states; X{t+1}];
    controls = [controls; U{t}];
    lams_l = [lams_l; YL{t}];
    lams_r = [lams_r; YR{t}];
    

    % Dynamic equations
    X_pred = X{t} + dt * [X{t}(2); ... %y
                         -g + (YL{t}+YR{t})/M; %vy
                          0]; %vtheta
    % Dynamic constraints.                
    h_x = [h_x; X_pred - X{t+1}];
    h_y_l{t} = [X_pred - X{t+1}];
    h_y_r{t} = [X_pred - X{t+1}];
    if (t>1)
        h_y_l{t-1} = [h_y_l{t-1}; X_pred-X{t+1}];
        h_y_r{t-1} = [h_y_r{t-1}; X_pred-X{t+1}];
    end


%     g_x = g_x;
    g_y_l{t} = []; % fill in later
    g_y_r{t} = []; % fill in later

    f_x = U{t}*U{t};
    
    f_y_l{t} = YL{t}'*YL{t};
    f_y_r{t} = YR{t}'*YR{t};
end

for t = 1:T

    if (t < T)
        p1 = X{t+2}(1) - l/2 * sin(X{t+2}(3));
        p2 = X{t+2}(1) + l/2 * sin(X{t+2}(3));
        g_y_l{t} = p1;
        g_y_r{t} = p2;
        vars_y_l{t} = [vars_y_l{t}; X{t+2}];
        vars_y_r{t} = [vars_y_r{t}; X{t+2}];
    else
        g_y_l{t} = [];
        g_y_r{t} = [];
    end
end

lam_x = MX.sym('lam_x', size(h_x,1));
mu_x = MX.sym('mu_x', size(g_x,1));
Lag_x = f_x; 
if size(h_x,1) > 0
    Lag_x = Lag_x - lam_x'*h_x;
end
if size(g_x,1) > 0
    Lag_x = Lag_x - mu_x'*g_x;
end
gLx = gradient(Lag_x, vars_x); % n + m
all_gLy_l = [];
all_gLy_r = [];
for t = 1:T
    lam_y_l{t} = MX.sym('lam_y_l', size(h_y_l{t},1));
    lam_y_r{t} = MX.sym('lam_y_r', size(h_y_r{t},1));
    mu_y_l{t} = MX.sym('mu_y_l', size(g_y_l{t},1));
    mu_y_r{t} = MX.sym('mu_y_r', size(g_y_r{t},1));

    Lag_y_l{t} = f_y_l{t};
    Lag_y_r{t} = f_y_r{t};

    if size(h_y_l{t},1) > 0
        Lag_y_l{t} = Lag_y_l{t} - lam_y_l{t}'*h_y_l{t};
    end
    if size(g_y_l{t},1) > 0
        Lag_y_l{t} = Lag_y_l{t} - mu_y_l{t}'*g_y_l{t};
    end
        if size(h_y_r{t},1) > 0
        Lag_y_r{t} = Lag_y_r{t} - lam_y_r{t}'*h_y_r{t};
    end
    if size(g_y_r{t},1) > 0
        Lag_y_r{t} = Lag_y_r{t} - mu_y_r{t}'*g_y_r{t};
    end
    all_gLy_l = [all_gLy_l; gradient(Lag_y_l{t}, vars_y_l{t})];
    all_gLy_r = [all_gLy_r; gradient(Lag_y_r{t}, vars_y_r{t})];
end


gL = [gLx; all_gLy_l; all_gLy_r];
% for t = 1:T
%     gL = [gL; gLy{t}; gLz{t}];
% end
gL = [gL; h_x; g_x];
all_g_y_l = [];
all_g_y_r = [];
for t = 1:T
    all_g_y_l = [all_g_y_l; g_y_l{t}];
    all_g_y_r = [all_g_y_r; g_y_r{t}];
%     gL = [gL; g_y{t}; g_z{t}];
end
gL = [gL; all_g_y_l; all_g_y_r];

all_vars = [states; controls; lams_l; lams_r; lam_x];
all_lam_y_l = [];
all_lam_y_r = [];
all_mu_y_l = [];
all_mu_y_r = [];
for t = 1:T
    all_lam_y_l = [all_lam_y_l; lam_y_l{t}];
    all_lam_y_r = [all_lam_y_r; lam_y_r{t}];
%     all_vars = [all_vars; lam_y{t}; lam_z{t}];
end
all_vars = [all_vars; all_lam_y_l; all_lam_y_r];
all_vars = [all_vars; mu_x];
for t = 1:T
    all_mu_y_l = [all_mu_y_l; mu_y_l{t}];
    all_mu_y_r = [all_mu_y_r; mu_y_r{t}];
%     all_vars = [all_vars; mu_y{t}; mu_z{t}];
end
all_vars = [all_vars; all_mu_y_l; all_mu_y_r];

HH = jacobian(gL, all_vars);
% spy(HH)

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});

extract_states = Function('extract_states', {all_vars}, {states});
extract_controls = Function('extract_controls', {all_vars}, {controls});
extract_lams_l = Function('extract_lams_l', {all_vars}, {lams_l});
extract_lams_r = Function('extract_lams_r', {all_vars}, {lams_r});

lo = [-inf_const*ones(size(states));
    -inf_const*ones(size(controls));
    -inf_const*ones(size(lams_l));
    -inf_const*ones(size(lams_r));
    -inf_const*ones(size(lam_x));
    -inf_const*ones(size(all_lam_y_l));
    -inf_const*ones(size(all_lam_y_r));
    zeros(size(mu_x));
    zeros(size(all_mu_y_l));
    zeros(size(all_mu_y_r))];
hi = [inf_const*ones(size(states));
    inf_const*ones(size(controls));
    inf_const*ones(size(lams_l));
    inf_const*ones(size(lams_r));
    inf_const*ones(size(lam_x));
    inf_const*ones(size(all_lam_y_l));
    inf_const*ones(size(all_lam_y_r));
    inf_const*ones(size(mu_x));
    inf_const*ones(size(all_mu_y_l));
    inf_const*ones(size(all_mu_y_r))];

z0 = -pi/4*ones(size(all_vars,1),1);

n = size(z0,1);
m = size([states; controls; lams_l; lams_r; lam_x; all_lam_y_l; all_lam_y_r], 1);
nnz = nnz(evalH(z0));
init_loc_a = 0;
init_loc_b = 0;
dim_state = size(X{1},1);
dim_control = size(U{1},1);
specs = [n;m;nnz;init_loc_a;init_loc_b;dim_state;dim_control;T;num_iterations];


writez(lo, '/Users/forrest/code/packages/pathlib-master/examples/C/lo.txt');
writez(hi, '/Users/forrest/code/packages/pathlib-master/examples/C/hi.txt');
writez(z0, '/Users/forrest/code/packages/pathlib-master/examples/C/z0.txt');
write_nums(specs, '/Users/forrest/code/packages/pathlib-master/examples/C/spec.txt');

opts = struct('with_header', true, 'casadi_int', 'int');

C = CodeGenerator('gen.c', opts);
C.add(evalG);
C.add(evalH);
C.generate();
disp 'Generated game!';

cmd = "cd /Users/forrest/code/packages/pathlib-master/examples/C && " + ...
      "export DYLD_LIBRARY_PATH=/Users/forrest/code/packages/pathlib-master/lib/osx && " + ...
      "export PATH_LICENSE_STRING='2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0' && " + ...
      "cp ~/Documents/MATLAB/fb_games/test_single_player/gen.c . && make && ./game_shared";
[status,cmdout] = system(cmd)

%%
z = readz('/Users/forrest/code/packages/pathlib-master/examples/C/z.csv');
states = full(extract_states(z));
controls = full(extract_controls(z));
lams_l = full(extract_lams_l(z));
lams_r = full(extract_lams_r(z));

y = states(1:3:end);
vy = states(2:3:end);
theta = states(3:3:end);

ry = sin(theta);
yr = y+l/2*ry;
yl = y-l/2*ry;

scale = 0.75;
figure; hold on;
plot([-5,5],[0,0],'-k')
axis([-scale,scale,-scale,scale]);
for t=1:T+1
    r = [cos(theta(t)); sin(theta(t))];
    p1 = [0; y(t)] - l/2 * r;
    p2 = [0; y(t)] + l/2 * r;
    p0{t} = plot([p1(1), p2(1)],[p1(2), p2(2)],'-k', 'linewidth', 4);
    axis([-scale,scale,-scale,scale]);
    pause(3*dt);
    delete(p0{t})
end


hold on;




