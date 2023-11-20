clear; clc; close all;
import casadi.*

% UNIFORM ROD FALLING!

inf_const = 999999999999.9;

l = 1;
M = 1;
g = 9.81;
I = 1/12 * M*l*l;
mu = 0.2;

dt = 0.1;
T = 10;
num_iterations = 1;

X{1} = MX.sym(['X_' num2str(1)], 12);

theta_offset = pi/2;

xa = [0; 4; 0; 0; 1; -3];
xb = [xa(1)+l/2+cos(theta_offset)*l/2; xa(2)+sin(theta_offset)*l/2; theta_offset; 0; 0; 0];

%     x   y   theta vx vy omega
x0 = [xa;xb];


% figure; hold on;
% ra = [cos(x0(3)); sin(x0(3))];
% p1a = [x0(1); x0(2)] - l/2 * ra;
% p2a = [x0(1); x0(2)] + l/2 * ra;
% rb = [cos(x0(9)); sin(x0(9))];
% p1b = [x0(7); x0(8)] - l/2 * rb;
% p2b = [x0(7); x0(8)] + l/2 * rb;
% prea = plot([p1a(1), p2a(1)],[p1a(2), p2a(2)],'-k', 'linewidth', 4);
% preb = plot([p1b(1), p2b(1)],[p1b(2), p2b(2)],'-k', 'linewidth', 4);
% axis([-7,7,-7,7]);
% axis('equal');

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
tqs = [];

for t = 1:T
    % Introduce symbolic vars.
    U{t} = MX.sym(['U_' num2str(t)], 1);
    % normal forces [l;r] (positive is up)
    YL{t} = MX.sym(['YL_' num2str(t)], 1);
    YR{t} = MX.sym(['YR_' num2str(t)], 1);
    % friction forces [l;r] (positive is right)
    ZL{t} = MX.sym(['Z_' num2str(t)], 1);
    ZR{t} = MX.sym(['Z_' num2str(t)], 1);
    TQ{t} = MX.sym(['TQ_' num2str(t)], 3);
    
    X{t+1} = MX.sym(['X_' num2str(t+1)], 12);

    vars_x = [vars_x; U{t}; X{t+1}];
    vars_y_l{t} = [YL{t}; X{t+1}];
    vars_y_r{t} = [YR{t}; X{t+1}];
    vars_z_l{t} = [ZL{t}; X{t+1}];
    vars_z_r{t} = [ZR{t}; X{t+1}];
    vars_tq{t} = [TQ{t}; X{t+1}];
   
    % Place vars in lists for easy access later.
    states = [states; X{t+1}];
    controls = [controls; U{t}];
    lams_l = [lams_l; YL{t}];
    lams_r = [lams_r; YR{t}];
    gams_l = [gams_l; ZL{t}];
    gams_r = [gams_r; ZR{t}];
    tqs = [tqs; TQ{t}];
    
    ra = [cos(X{t}(3));
        sin(X{t}(3))];
    rb = [cos(X{t}(9));
        sin(X{t}(9))];
    
    if (t==1)
        push = [-3; -3];
    else
        push = [0;0];
    end
    
    % Dynamic equations
    X_pred = X{t} + dt * [X{t}(4); ... %x
                          X{t}(5); ... %y
                          X{t}(6); ... %t
                          (push(1)+TQ{t}(1))/M; ... %vx
                          (push(2)+TQ{t}(2))/M; ... %vy
                          l/(2*I)*(push(1)*ra(2)-push(2)*ra(1)-TQ{t}(1)*ra(2)+TQ{t}(2)*ra(1));
                          X{t}(10); ... %x
                          X{t}(11); ... %y
                          X{t}(12); ... %t
                          -TQ{t}(1)/M; ... %vx
                          -TQ{t}(2)/M; ... %vy
                          l/(2*I)*(-TQ{t}(1)*rb(2)+TQ{t}(2)*rb(1))];
    X_pred(5) = X_pred(5) + dt * (-g);
    X_pred(11) = X_pred(11) + dt*(-g);
        
    
    % Dynamic constraints.                
    h_x = [h_x; X_pred - X{t+1}];
    h_tq{t} = [X_pred - X{t+1}];
    h_y_l{t} = [X_pred - X{t+1}];
    h_y_r{t} = [X_pred - X{t+1}];
    h_z_l{t} = [X_pred - X{t+1}];
    h_z_r{t} = [X_pred - X{t+1}];
    if (t>1)
        h_tq{t-1} = [h_tq{t-1}; X_pred-X{t+1}];
        h_y_l{t-1} = [h_y_l{t-1}; X_pred-X{t+1}];
        h_y_r{t-1} = [h_y_r{t-1}; X_pred-X{t+1}];
    end


%     g_x = g_x;
    g_y_l{t} = []; % fill in later
    g_y_r{t} = []; % fill in later
    g_z_l{t} = [];%0.00000001+YL{t}*YL{t}*mu^2 - ZL{t}*ZL{t};
    g_z_r{t} = [];%0.00000001+YR{t}*YR{t}*mu^2 - ZR{t}*ZR{t};

    f_x = U{t}*U{t};% + X{t+1}(3)*X{t+1}(3);
    
    f_y_l{t} = YL{t}'*YL{t};
    f_y_r{t} = YR{t}'*YR{t};
    
    % TODO! This should be velocity of left and right points.
%     f_z_l{t} = (X{t+1}(4) + l/2 * sin(X{t+1}(3))*X{t+1}(6))^2;
%     f_z_r{t} = (X{t+1}(4) - l/2 * sin(X{t+1}(3))*X{t+1}(6))^2;
    f_z_l{t} = ZL{t}*ZL{t};
    f_z_r{t} = ZR{t}*ZR{t};
    
    f_tq{t} = TQ{t}'*TQ{t};

end

for t = 1:T

    if (t < T)
        r = [cos(X{t+2}(3)); sin(X{t+2}(3))];
        p1 = [X{t+2}(1);X{t+2}(2)] - l/2 * r;
        p2 = [X{t+2}(1);X{t+2}(2)] + l/2 * r;
%         g_y_l{t} = p1(2);
%         g_y_r{t} = p2(2);
        vars_y_l{t} = [vars_y_l{t}; X{t+2}];
        vars_y_r{t} = [vars_y_r{t}; X{t+2}];
        vars_tq{t} = [vars_tq{t}; X{t+2}];
        hh_tq{t} = [X{t+2}(1)+l/2*cos(X{t+2}(3)) - (X{t+2}(7)-l/2*cos(X{t+2}(9)));
                    X{t+2}(2)+l/2*sin(X{t+2}(3)) - (X{t+2}(8)-l/2*sin(X{t+2}(9)))];
                    %X{t+2}(3) - X{t+2}(9)+theta_offset];
    else
        g_y_l{t} = [];
        g_y_r{t} = [];
        hh_tq{t} = [];
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
all_gLtq = [];
for t = 1:T
    lam_y_l{t} = MX.sym('lam_y_l', size(h_y_l{t},1));
    lam_y_r{t} = MX.sym('lam_y_r', size(h_y_r{t},1));
    lam_z_l{t} = MX.sym('lam_z_l', size(h_z_l{t},1));
    lam_z_r{t} = MX.sym('lam_z_r', size(h_z_r{t},1));
    lam_tq{t} = MX.sym('lam_tq', size([h_tq{t}; hh_tq{t}],1));
    mu_y_l{t} = MX.sym('mu_y_l', size(g_y_l{t},1));
    mu_y_r{t} = MX.sym('mu_y_r', size(g_y_r{t},1));
    mu_z_l{t} = MX.sym('mu_z_l', size(g_z_l{t},1));
    mu_z_r{t} = MX.sym('mu_z_r', size(g_z_r{t},1));

    Lag_y_l{t} = f_y_l{t};
    Lag_y_r{t} = f_y_r{t};
    Lag_z_l{t} = f_z_l{t};
    Lag_z_r{t} = f_z_r{t};
    Lag_tq{t} = f_tq{t};


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
    if size([h_tq{t}; hh_tq{t}],1) > 0
        Lag_tq{t} = Lag_tq{t} - lam_tq{t}'*[h_tq{t}; hh_tq{t}];
    end
    all_gLy_l = [all_gLy_l; gradient(Lag_y_l{t}, vars_y_l{t})];
    all_gLy_r = [all_gLy_r; gradient(Lag_y_r{t}, vars_y_r{t})];
    all_gLz_l = [all_gLz_l; gradient(Lag_z_l{t}, vars_z_l{t})];
    all_gLz_r = [all_gLz_r; gradient(Lag_z_r{t}, vars_z_r{t})];
    all_gLtq = [all_gLtq; gradient(Lag_tq{t}, vars_tq{t})];
%     gLy{t} = gradient(Lag_y{t}, vars_y{t}); % n + l
%     gLz{t} = gradient(Lag_z{t}, vars_z{t}); % n + g
end


gL = [gLx; all_gLy_l; all_gLy_r; all_gLz_l; all_gLz_r; all_gLtq];
% for t = 1:T
%     gL = [gL; gLy{t}; gLz{t}];
% end
gL = [gL; h_x; g_x];
all_g_y_l = [];
all_g_y_r = [];
all_g_z_l = [];
all_g_z_r = [];
all_hh_tq = [];
for t = 1:T
    all_g_y_l = [all_g_y_l; g_y_l{t}];
    all_g_y_r = [all_g_y_r; g_y_r{t}];
    all_g_z_l = [all_g_z_l; g_z_l{t}];
    all_g_z_r = [all_g_z_r; g_z_r{t}];
    all_hh_tq = [all_hh_tq; hh_tq{t}];
%     gL = [gL; g_y{t}; g_z{t}];
end
gL = [gL; all_hh_tq; all_g_y_l; all_g_y_r; all_g_z_l; all_g_z_r];

all_vars = [states; controls; lams_l; lams_r; gams_l; gams_r; tqs; lam_x];
all_lam_y_l = [];
all_lam_y_r = [];
all_lam_z_l = [];
all_lam_z_r = [];
all_lam_tq = [];
all_mu_y_l = [];
all_mu_y_r = [];
all_mu_z_l = [];
all_mu_z_r = [];
for t = 1:T
    all_lam_y_l = [all_lam_y_l; lam_y_l{t}];
    all_lam_y_r = [all_lam_y_r; lam_y_r{t}];
    all_lam_z_l = [all_lam_z_l; lam_z_l{t}];
    all_lam_z_r = [all_lam_z_r; lam_z_r{t}];
    all_lam_tq = [all_lam_tq; lam_tq{t}];
%     all_vars = [all_vars; lam_y{t}; lam_z{t}];
end
all_vars = [all_vars; all_lam_y_l; all_lam_y_r; all_lam_z_l; all_lam_z_r; all_lam_tq];
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
extract_tqs = Function('extract_tqs', {all_vars}, {tqs});


lo = [-inf_const*ones(size(states));
    -inf_const*ones(size(controls));
    -inf_const*ones(size(lams_l));
    -inf_const*ones(size(lams_r));
    -inf_const*ones(size(gams_l));
    -inf_const*ones(size(gams_r));
    -inf_const*ones(size(tqs));
    -inf_const*ones(size(lam_x));
    -inf_const*ones(size(all_lam_y_l));
    -inf_const*ones(size(all_lam_y_r));
    -inf_const*ones(size(all_lam_z_l));
    -inf_const*ones(size(all_lam_z_r));
    -inf_const*ones(size(all_lam_tq));
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
    inf_const*ones(size(tqs));
    inf_const*ones(size(lam_x));
    inf_const*ones(size(all_lam_y_l));
    inf_const*ones(size(all_lam_y_r));
    inf_const*ones(size(all_lam_z_l));
    inf_const*ones(size(all_lam_z_r));
    inf_const*ones(size(all_lam_tq));
    inf_const*ones(size(mu_x));
    inf_const*ones(size(all_mu_y_l));
    inf_const*ones(size(all_mu_y_r));
    inf_const*ones(size(all_mu_z_l));
    inf_const*ones(size(all_mu_z_r))];

z0 = rand(size(all_vars,1),1);
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

xa = states(1:12:end);
ya = states(2:12:end);
thetaa = states(3:12:end);
vxa = states(4:12:end);
vya = states(5:12:end);
omegaa = states(6:12:end);
xb = states(7:12:end);
yb = states(8:12:end);
thetab = states(9:12:end);
vxb = states(10:12:end);
vyb = states(11:12:end);
omegab = states(12:12:end);

rxa = cos(thetaa);
rya = sin(thetaa);
yra = ya+l/2*rya;
yla = ya-l/2*rya;

rxb = cos(thetab);
ryb = sin(thetab);
yrb = yb+l/2*ryb;
ylb = yb-l/2*ryb;


scale = 10;
figure; hold on;
plot([min(xa)-5,max(xa)+5],[0,0],'-k')
axis([-scale,scale,-scale,scale]);
for t=1:T+1
    ra = [cos(thetaa(t)); sin(thetaa(t))];
    p1a = [xa(t); ya(t)] - l/2 * ra;
    p2a = [xa(t); ya(t)] + l/2 * ra;
    rb = [cos(thetab(t)); sin(thetab(t))];
    p1b = [xb(t); yb(t)] - l/2 * rb;
    p2b = [xb(t); yb(t)] + l/2 * rb;
    p0a{t} = plot([p1a(1), p2a(1)],[p1a(2), p2a(2)],'-g', 'linewidth', 4);
    p0b{t} = plot([p1b(1), p2b(1)],[p1b(2), p2b(2)],'-b', 'linewidth', 4);
    axis([-scale,scale,-scale,scale]);
    pause(4*dt);
    delete(p0a{t})
    delete(p0b{t})
end


hold on;




