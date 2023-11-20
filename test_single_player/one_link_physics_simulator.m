clear; clc; close all;
import casadi.*

% UNIFORM ROD FALLING!

inf_const = 999999999999.9;

l = 1;
M = 1;
g = 9.81;
I = 1/12 * M*l*l;
mu = 0.1;

dt = 0.05;
T = 10;
num_iterations = 1;

X{1} = MX.sym(['X_' num2str(1)], 6);

t0 = -pi/3;
%     x   y   theta vx vy omega
x0 = [0; 1.75; t0; 3; -5; 0];

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
    % normal forces [l;r] (positive is up)
    YL{t} = MX.sym(['YL_' num2str(t)], 1);
    YR{t} = MX.sym(['YR_' num2str(t)], 1);
    % friction forces [l;r] (positive is right)
    ZL{t} = MX.sym(['Z_' num2str(t)], 1);
    ZR{t} = MX.sym(['Z_' num2str(t)], 1);
    
    X{t+1} = MX.sym(['X_' num2str(t+1)], 6);

    vars_x = [vars_x; U{t}; X{t+1}];
%     vars_y_l{t} = YL{t}; 
%     vars_y_r{t} = YR{t}; 
    vars_y_l{t} = [YL{t}; X{t+1}];
    vars_y_r{t} = [YR{t}; X{t+1}];
    vars_z_l{t} = [ZL{t}; X{t+1}];
    vars_z_r{t} = [ZR{t}; X{t+1}];
   
    % Place vars in lists for easy access later.
    states = [states; X{t+1}];
    controls = [controls; U{t}];
    lams_l = [lams_l; YL{t}];
    lams_r = [lams_r; YR{t}];
    gams_l = [gams_l; ZL{t}];
    gams_r = [gams_r; ZR{t}];

    r = [cos(X{t}(3));
        sin(X{t}(3))];

    % Dynamic equations
    if (t < T) 
        X_pred = X{t} + dt * [X{t}(4); ... %x
                              X{t}(5); ... %y
                              X{t}(6); ... %t
                              (ZL{t}+ZR{t})/M; ... %vx
                              (YL{t}+YR{t})/M; ... %vy
                              l/(2*I)*(ZL{t}*r(2)-YL{t}*r(1) - ZR{t}*r(2) +YR{t}*r(1))];
    else
        X_pred = X{t} + dt * [X{t}(4); ... %x
                              X{t}(5); ... %y
                              X{t}(6); ... %t
                              (ZL{t}+ZR{t})/M; ... %vx
                              (YL{t}+YR{t})/M; ... %vy
                              l/(2*I)*(ZL{t}*r(2)-YL{t}*r(1) - ZR{t}*r(2) +YR{t}*r(1))];
    end
    X_pred(5) = X_pred(5) + dt * (-g);
    
    % Dynamic constraints.                
    h_x = [h_x; X_pred - X{t+1}];
    h_y_l{t} = [X_pred - X{t+1}];
    h_y_r{t} = [X_pred - X{t+1}];

    h_z_l{t} = [X_pred - X{t+1}];
    h_z_r{t} = [X_pred - X{t+1}];
    if (t>1)
        h_y_l{t-1} = [h_y_l{t-1}; X_pred-X{t+1}];
        h_y_r{t-1} = [h_y_r{t-1}; X_pred-X{t+1}];
    end


%     g_x = g_x;
    g_y_l{t} = []; % fill in later
    g_y_r{t} = []; % fill in later
    g_z_l{t} = 0.0000001+YL{t}*YL{t}*mu^2 - ZL{t}*ZL{t};
    g_z_r{t} = 0.0000001+YR{t}*YR{t}*mu^2 - ZR{t}*ZR{t};
%     g_z_l{t} = [];
%     g_z_r{t} = [];

    f_x = U{t}*U{t};% + X{t+1}(3)*X{t+1}(3);
    
    f_y_l{t} = YL{t}'*YL{t};
    f_y_r{t} = YR{t}'*YR{t};
    
    % TODO! This should be velocity of left and right points.
    f_z_l{t} = (X{t+1}(4) + l/2 * sin(X{t+1}(3))*X{t+1}(6))^2;
    f_z_r{t} = (X{t+1}(4) - l/2 * sin(X{t+1}(3))*X{t+1}(6))^2;
%     f_z_l{t} = ZL{t}*ZL{t};
%     f_z_r{t} = ZR{t}*ZR{t};

end

for t = 1:T

    if (t < T)
        r = [cos(X{t+2}(3)); sin(X{t+2}(3))];
        pc = [X{t+2}(1);X{t+2}(2)];
        pr = pc + l/2*r;
        pl = pc - l/2*r;

        g_y_l{t} = pl(2);
        g_y_r{t} = pr(2);
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
all_gLz_l = [];
all_gLz_r = [];
for t = 1:T
    lam_y_l{t} = MX.sym('lam_y_l', size(h_y_l{t},1));
    lam_y_r{t} = MX.sym('lam_y_r', size(h_y_r{t},1));
    lam_z_l{t} = MX.sym('lam_z_l', size(h_z_l{t},1));
    lam_z_r{t} = MX.sym('lam_z_r', size(h_z_r{t},1));
    mu_y_l{t} = MX.sym('mu_y_l', size(g_y_l{t},1));
    mu_y_r{t} = MX.sym('mu_y_r', size(g_y_r{t},1));
    mu_z_l{t} = MX.sym('mu_z_l', size(g_z_l{t},1));
    mu_z_r{t} = MX.sym('mu_z_r', size(g_z_r{t},1));

    Lag_y_l{t} = f_y_l{t};
    Lag_y_r{t} = f_y_r{t};
    Lag_z_l{t} = f_z_l{t};
    Lag_z_r{t} = f_z_r{t};


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
    if size(h_z_l{t},1) > 0
        Lag_z_l{t} = Lag_z_l{t} - lam_z_l{t}'*h_z_l{t};
    end
    if size(g_z_l{t},1) > 0
        Lag_z_l{t} = Lag_z_l{t} - mu_z_l{t}'*g_z_l{t};
    end
    if size(h_z_r{t},1) > 0
        Lag_z_r{t} = Lag_z_r{t} - lam_z_r{t}'*h_z_r{t};
    end
    if size(g_z_r{t},1) > 0
        Lag_z_r{t} = Lag_z_r{t} - mu_z_r{t}'*g_z_r{t};
    end
    all_gLy_l = [all_gLy_l; gradient(Lag_y_l{t}, vars_y_l{t})];
    all_gLy_r = [all_gLy_r; gradient(Lag_y_r{t}, vars_y_r{t})];
    all_gLz_l = [all_gLz_l; gradient(Lag_z_l{t}, vars_z_l{t})];
    all_gLz_r = [all_gLz_r; gradient(Lag_z_r{t}, vars_z_r{t})];
%     gLy{t} = gradient(Lag_y{t}, vars_y{t}); % n + l
%     gLz{t} = gradient(Lag_z{t}, vars_z{t}); % n + g
end


gL = [gLx; all_gLy_l; all_gLy_r; all_gLz_l; all_gLz_r];
% for t = 1:T
%     gL = [gL; gLy{t}; gLz{t}];
% end
gL = [gL; h_x; g_x];
all_g_y_l = [];
all_g_y_r = [];
all_g_z_l = [];
all_g_z_r = [];
for t = 1:T
    all_g_y_l = [all_g_y_l; g_y_l{t}];
    all_g_y_r = [all_g_y_r; g_y_r{t}];
    all_g_z_l = [all_g_z_l; g_z_l{t}];
    all_g_z_r = [all_g_z_r; g_z_r{t}];
%     gL = [gL; g_y{t}; g_z{t}];
end
gL = [gL; all_g_y_l; all_g_y_r; all_g_z_l; all_g_z_r];

all_vars = [states; controls; lams_l; lams_r; gams_l; gams_r; lam_x];
all_lam_y_l = [];
all_lam_y_r = [];
all_lam_z_l = [];
all_lam_z_r = [];
all_mu_y_l = [];
all_mu_y_r = [];
all_mu_z_l = [];
all_mu_z_r = [];
for t = 1:T
    all_lam_y_l = [all_lam_y_l; lam_y_l{t}];
    all_lam_y_r = [all_lam_y_r; lam_y_r{t}];
    all_lam_z_l = [all_lam_z_l; lam_z_l{t}];
    all_lam_z_r = [all_lam_z_r; lam_z_r{t}];
%     all_vars = [all_vars; lam_y{t}; lam_z{t}];
end
all_vars = [all_vars; all_lam_y_l; all_lam_y_r; all_lam_z_l; all_lam_z_r];
all_vars = [all_vars; mu_x];
for t = 1:T
    all_mu_y_l = [all_mu_y_l; mu_y_l{t}];
    all_mu_y_r = [all_mu_y_r; mu_y_r{t}];
    all_mu_z_l = [all_mu_z_l; mu_z_l{t}];
    all_mu_z_r = [all_mu_z_r; mu_z_r{t}];
%     all_vars = [all_vars; mu_y{t}; mu_z{t}];
end
all_vars = [all_vars; all_mu_y_l; all_mu_y_r; all_mu_z_l; all_mu_z_r];

HH = jacobian(gL, all_vars);
% spy(HH)

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});

extract_states = Function('extract_states', {all_vars}, {states});
extract_controls = Function('extract_controls', {all_vars}, {controls});
extract_lams_l = Function('extract_lams_l', {all_vars}, {lams_l});
extract_lams_r = Function('extract_lams_r', {all_vars}, {lams_r});
extract_gams_l = Function('extract_gams_l', {all_vars}, {gams_l});
extract_gams_r = Function('extract_gams_r', {all_vars}, {gams_r});


lo = [-inf_const*ones(size(states));
    -inf_const*ones(size(controls));
    -inf_const*ones(size(lams_l));
    -inf_const*ones(size(lams_r));
    -inf_const*ones(size(gams_l));
    -inf_const*ones(size(gams_r));
    -inf_const*ones(size(lam_x));
    -inf_const*ones(size(all_lam_y_l));
    -inf_const*ones(size(all_lam_y_r));
    -inf_const*ones(size(all_lam_z_l));
    -inf_const*ones(size(all_lam_z_r));
    zeros(size(mu_x));
    zeros(size(all_mu_y_l));
    zeros(size(all_mu_y_r));
    zeros(size(all_mu_z_l));
    zeros(size(all_mu_z_r))];
hi = [inf_const*ones(size(states));
    inf_const*ones(size(controls));
    inf_const*ones(size(lams_l));
    inf_const*ones(size(lams_r));
    inf_const*ones(size(gams_l));
    inf_const*ones(size(gams_r));
    inf_const*ones(size(lam_x));
    inf_const*ones(size(all_lam_y_l));
    inf_const*ones(size(all_lam_y_r));
    inf_const*ones(size(all_lam_z_l));
    inf_const*ones(size(all_lam_z_r));
    inf_const*ones(size(mu_x));
    inf_const*ones(size(all_mu_y_l));
    inf_const*ones(size(all_mu_y_r));
    inf_const*ones(size(all_mu_z_l));
    inf_const*ones(size(all_mu_z_r))];

cur = x0;
states0 = x0;
for t = 1:T
    cur = cur + dt* [cur(4); ... %x
                          cur(5); ... %y
                          cur(6); ... %t
                          0; ... %vx
                          -g; ... %vy
                          0];
    if (cur(2) < 0) 
        cur(2) = 0;
        cur(5) = 0;
    end
    states0 = [states0; cur];
end


z0 = rand(size(all_vars,1),1);
% z0(1:size(states,1)) = states0;
% z0 = load('good_soln_z.mat');
% z0 = z0.z;
% z0 = full(z);

n = size(z0,1);
m = size([states; controls; lams_l; lams_r; gams_l; gams_r; lam_x; all_lam_y_l; all_lam_y_r; all_lam_z_l; all_lam_z_r], 1);
nnz = nnz(evalH(z0));
init_loc_a = 0;
init_loc_b = 0;
dim_state = size(X{1},1);
dim_control = size(U{1},1);
specs = [n;m;nnz;init_loc_a;init_loc_b;dim_state;dim_control;T;num_iterations;dt];


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
gams_l = full(extract_gams_l(z));
gams_r = full(extract_gams_r(z));

x = states(1:6:end);
y = states(2:6:end);
theta = states(3:6:end);
vx = states(4:6:end);
vy = states(5:6:end);
omega = states(6:6:end);

rx = cos(theta);
ry = sin(theta);
yr = y(3:end)+l/2*ry(1:end-2);
yl = y(3:end)-l/2*ry(1:end-2);


scale = 4;
figure; hold on;
plot([min(x)-5,max(x)+5],[0,0],'-k')
axis([-scale,scale,-scale,scale]);
for t=1:T+1
    r = [cos(theta(t)); sin(theta(t))];
    p1 = [x(t); y(t)] - l/2 * r;
    p2 = [x(t); y(t)] + l/2 * r;
    p0{t} = plot([p1(1), p2(1)],[p1(2), p2(2)],'-k', 'linewidth', 4);
    axis([-scale,scale,-scale,scale]);
    pause(1*dt);
    delete(p0{t})
end


hold on;




