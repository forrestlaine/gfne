clear; clc; close all;
import casadi.*

inf_const = 999999999999.9;

dt = 0.25;
T = 20;
num_iterations = 80;

X_a{1} = MX.sym(['X_a_' num2str(0)], 4);
X_b{1} = MX.sym(['X_b_' num2str(0)], 4);

a_bound_a = inf_const; % m/s^2
a_bound_b = inf_const;
v_bound = inf_const;

x0_a = [10;-8;0;0];
x0_b = [3;-4;0;0];

f_a = 0;
f_b = 0;

f2_a = 0;
f2_b = 0;

h_a = [];
h_b = [];
h2_a = [];
h2_b = [];

lb_a = x0_a;
lb_b = x0_b;

ub_a = x0_a;
ub_b = x0_b;

g_a = [];
g_b = [];
g2_a = [];
g2_b = [];

vars_a = X_a{1};
vars_b = X_b{1};

vars2_a = [];
vars2_b = [];

p_a = X_a{1};
p_b = X_b{1};

q_a = [];
q_b = [];

for t = 1:T
    % Introduce symbolic vars.
    U_a{t} = MX.sym(['U_a_' num2str(t)], 2);
    U_b{t} = MX.sym(['U_b_' num2str(t)], 2);

    X_a{t+1} = MX.sym(['X_a_' num2str(t)], 4);
    X_b{t+1} = MX.sym(['X_b_' num2str(t)], 4);
   

    % Place vars in lists for easy access later.
    p_a = [p_a; X_a{t+1}];
    p_b = [p_b; X_b{t+1}];
    
    q_a = [q_a; U_a{t}];
    q_b = [q_b; U_b{t}];
    
    % Dynamic equations
    X_a_pred = X_a{t} + dt * [X_a{t}(3); ...
                              X_a{t}(4); ...
                              U_a{t}(1); ...
                              U_a{t}(2)];
    X_b_pred = X_b{t} + dt * [X_b{t}(3); ...
                              X_b{t}(4); ...
                              U_b{t}(1); ...
                              U_b{t}(2)];             

    % Primal variable bounds.
    lb_a = [lb_a; -a_bound_a; -a_bound_a; -inf_const; -inf_const; -v_bound; -v_bound];
    ub_a = [ub_a;  a_bound_a; a_bound_a;  inf_const; inf_const; v_bound; v_bound];
    lb_b = [lb_b; -a_bound_b; -a_bound_b; -inf_const; -inf_const; -v_bound; -v_bound];
    ub_b = [ub_b;  a_bound_b; a_bound_b; inf_const; inf_const; v_bound; v_bound];
    
    eqs_a = [X_a_pred - X_a{t+1}];
    eqs_b = [X_b_pred - X_b{t+1}];
    
    ineqs_a = [(X_a{t+1}(1:2))'*(X_a{t+1}(1:2)) - 9;
        (X_a{t+1}(1:2)-X_b{t+1}(1:2))'*(X_a{t+1}(1:2)-X_b{t+1}(1:2)) - 4];
    ineqs_b = [(X_b{t+1}(1:2))'*(X_b{t+1}(1:2)) - 9];
    
    cost_a = 1*(U_a{t}(1:2)'*U_a{t}(1:2)) + ... 
        1*(X_a{t+1}(1:2)'*X_a{t+1}(1:2));
    cost_b = 1*(U_b{t}(1:2)'*U_b{t}(1:2)) + ...
         1*(X_b{t+1}(1:2) - X_a{t+1}(1:2))' * (X_b{t+1}(1:2) - X_a{t+1}(1:2));
    
    
    vars_a = [vars_a; U_a{t}; X_a{t+1}];
    vars_b = [vars_b; U_b{t}; X_b{t+1}];
    
    if (t <= T/2)    
        f_a = f_a + cost_a;
        f_b = f_b + cost_b;

        h_a = [h_a; eqs_a];
        h_b = [h_b; eqs_b];

        g_a = [g_a; 
            ineqs_a];
        g_b = [g_b; 
            ineqs_b];
    end
     
    if (t > T/2)
        vars2_a = [vars2_a; U_a{t}; X_a{t+1}];
        vars2_b = [vars2_b; U_b{t}; X_b{t+1}];
        f2_a = f2_a + cost_a;
        f2_b = f2_b + cost_b;

        h2_a = [h2_a; eqs_a];
        h2_b = [h2_b; eqs_b];

        g2_a = [g2_a; 
            ineqs_a];
        g2_b = [g2_b; 
            ineqs_b];
    end
end

lam2_a = MX.sym('lam2_a', size(h2_a,1));
lam2_b = MX.sym('lam2_b', size(h2_b,1));
mu2_a = MX.sym('mu2_a', size(g2_a,1));
mu2_b = MX.sym('mu2_b', size(g2_b,1));

Lag2_a = f2_a;
Lag2_b = f2_b;
if size(h2_a,1) > 0
    Lag2_a = Lag2_a - lam2_a'*h2_a;
end
if size(g2_a,1) > 0
    Lag2_a = Lag2_a - mu2_a'*g2_a;
end
if size(h2_b,1) > 0
    Lag2_b = Lag2_b - lam2_b'*h2_b;
end
if size(g2_b,1) > 0
    Lag2_b = Lag2_b - mu2_b'*g2_b;
end

gL2a = gradient(Lag2_a, vars2_a);
gL2b = gradient(Lag2_b, vars2_b);

shared_vars = [vars2_a; vars2_b; lam2_a; lam2_b; mu2_a; mu2_b];
shared_eq = [gL2a; gL2b; h2_a; h2_b];
shared_ineq = [g2_a; g2_b];

h_a = [h_a; shared_eq];
h_b = [h_b; shared_eq];
g_a = [g_a; shared_ineq];
g_b = [g_b; shared_ineq];

lam_a = MX.sym('lam_a', size(h_a,1));
lam_b = MX.sym('lam_b', size(h_b,1));

mu_a = MX.sym('mu_a', size(g_a,1));
mu_b = MX.sym('mu_b', size(g_b,1));

Lag_a = f_a; 
Lag_b = f_b;

if size(h_a,1) > 0
    Lag_a = Lag_a - lam_a'*h_a;
end
if size(g_a,1) > 0
    Lag_a = Lag_a - mu_a'*g_a;
end
if size(h_b,1) > 0
    Lag_b = Lag_b - lam_b'*h_b;
end
if size(g_b,1) > 0
    Lag_b = Lag_b - mu_b'*g_b;
end

vars_a = [vars_a; slacks_a];
vars_b = [vars_b; slacks_b];

dvars_a = [vars_a; shared_vars];
dvars_b = [vars_b; shared_vars];

gLa = gradient(Lag_a, dvars_a);
gLb = gradient(Lag_b, dvars_b);

gL = [gLa; gLb; h_a; h_b; shared_eq; g_a; g_b];

all_vars = [vars_a; vars_b; lam_a; lam_b; shared_vars; mu_a; mu_b];
HH = jacobian(gL, all_vars);

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});
evalFa = Function('obj_a', {all_vars}, {f_a});
evalFb = Function('obj_b', {all_vars}, {f_b});
evalha = Function('const_a', {all_vars}, {h_a});
evalhb = Function('const_b', {all_vars}, {h_b});

evalga = Function('const_a', {all_vars}, {g_a});
evalgb = Function('const_b', {all_vars}, {g_b});

extract_pa = Function('extract_pa', {all_vars}, {p_a});
extract_pb = Function('extract_pb', {all_vars}, {p_b});

extract_qa = Function('extract_qa', {all_vars}, {q_a});
extract_qb = Function('extract_qb', {all_vars}, {q_b});

extract_lam_a = Function('extract_mu_a', {all_vars}, {lam_a});
extract_lam_b = Function('extract_mu_b', {all_vars}, {lam_b});

extract_mu_a = Function('extract_mu_a', {all_vars}, {mu_a});
extract_mu_b = Function('extract_mu_b', {all_vars}, {mu_b});


ones_a = ones(size(vars_a,1),1);
ones_lam_a = ones(size(lam_a,1),1);
ones_mu_a = ones(size(mu_a,1),1);

ones_b = ones(size(vars_b,1),1);
ones_lam_b = ones(size(lam_b,1),1);
ones_mu_b = ones(size(mu_b,1),1);

lo = [lb_a; lb_b; 
    -inf_const*ones_lam_a;
    -inf_const*ones_lam_b;
    0*ones_mu_a;
    0*ones_mu_b];
hi = [ub_a; ub_b; 
    inf_const*ones_lam_a;
    inf_const*ones_lam_b;
    inf_const*ones_mu_a;
    inf_const*ones_mu_b];

z0 = zeros(size(all_vars,1),1);
xx_a = x0_a;
xx_b = x0_b;
init_a = x0_a;
init_b = x0_b;

for t = 1:T
    u_a = [0;0];
    u_b = [0;0];
    xx_a = xx_a + dt * [xx_a(3);xx_a(4);u_a];
    init_a = [init_a; u_a; xx_a];
    xx_b = xx_b + dt * [xx_b(3);xx_b(4);u_b];
    init_b = [init_b; u_b; xx_b];
end

XA = init_a(1:6:end);
YA = init_a(2:6:end);

XB = init_b(1:6:end);
YB = init_b(2:6:end);

S = size(vars_a,1);
z0(1:S) = init_a;
z0(S+1:2*S) = init_b;

n = size(z0,1);
m = size([vars_a; vars_b;lam_a; lam_b],1);
nnz = nnz(evalH(z0));
init_loc_a = 0;
init_loc_b = size(vars_a,1);
dim_state = size(X_a{1},1);
dim_control = size(U_a{1},1);
specs = [n;m;nnz;init_loc_a;init_loc_b;dim_state;dim_control;T;num_iterations; dt];


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
ZZ = readz('/Users/forrest/code/packages/pathlib-master/examples/C/z.csv');

figure; hold on;
viscircles([0,0],2);
for iter = 1:num_iterations
    z = ZZ((iter-1)*n+1:iter*n);

    PA = full(extract_pa(z));
    PB = full(extract_pb(z));

    QA = full(extract_qa(z));
    QB = full(extract_qb(z));

    XA = PA(1:4:end);
    YA = PA(2:4:end);
    AA = QA(2:2:end);
    LA = QA(1:2:end);

    XB = PB(1:4:end);
    YB = PB(2:4:end);
    AB = QB(2:2:end);
    LB = QB(1:2:end);
   
    path_a = plot(XA(:,1), YA(:,1),'or', 'MarkerSize', 1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
    path_b = plot(XB(:,1),YB(:,1),'or', 'MarkerSize', 1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'r');

%     spot_a = plot(XA(1,1), YA(1,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     spot_b = plot(XB(1,1),YB(1,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    cur_a = viscircles([XA(1,1),YA(1,1)],1,'Color', 'b');
    cur_b = viscircles([XB(1,1),YB(1,1)],1,'Color', 'r');
    time = text(XA(1,1)+15, YA(1,1)+15,string(iter*dt));
    axis([XA(1,1)-20 XA(1,1)+20 YA(1,1)-20 YA(1,1)+20]);
    pause(dt);
    delete(path_a);
    delete(path_b);
    delete(cur_a);
    delete(cur_b);
    delete(time);
end


