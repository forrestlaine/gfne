clear; clc; close all;
import casadi.*

inf_const = 999999999999.9;

dt = 0.1;
T = 10;

X_a{1} = MX.sym(['X_a_' num2str(0)], 3);
X_b{1} = MX.sym(['X_b_' num2str(0)], 3);
X_c{1} = MX.sym(['X_c_' num2str(0)], 3);

% Referee primal vars.
P_ab{1} = MX.sym(['P_ab_' num2str(0)], 2);
P_ba{1} = MX.sym(['P_ba_' num2str(0)], 2);
P_ac{1} = MX.sym(['P_ac_' num2str(0)], 2);
P_ca{1} = MX.sym(['P_ca_' num2str(0)], 2);
P_bc{1} = MX.sym(['P_bc_' num2str(0)], 2);
P_cb{1} = MX.sym(['P_cb_' num2str(0)], 2);

VBASE = 30;

va0 = 30-VBASE;
vb0 = 30-VBASE;
vc0 = 30-VBASE;

va = 33-VBASE;
vb = 31-VBASE;
vc = 29-VBASE;

half_w = 1.5;
half_h = 5.0;
shear = 0.2;

lat_bound = inf_const; % m/s
a_bound = inf_const; % m/s^2

x0_a = [-10;0;va0];
x0_b = [10;0;va0];
x0_c = [30;0;va0];

f_a = 0;
f_b = 0;
f_c = 0;
f_r = 0;
f_r = 1*(P_ab{1}-P_ba{1})'*(P_ab{1}-P_ba{1}) + ...
      1*(P_ac{1}-P_ca{1})'*(P_ac{1}-P_ca{1}) + ...
      1*(P_bc{1}-P_cb{1})'*(P_bc{1}-P_cb{1});

h_a = [];
h_b = [];
h_c = [];
h_r = [];
h_shared = [];

lb_a = x0_a;
lb_b = x0_b;
lb_c = x0_c;
ub_a = x0_a;
ub_b = x0_b;
ub_c = x0_c;

g_a = [];
g_b = [];
g_c = [];
g_r = [];
g_r = [1 - ((P_ab{1}(2)-X_a{1}(2))/half_w)^4 - (((P_ab{1}(1)-X_a{1}(1)) + shear*(P_ab{1}(2)-X_a{1}(2)))/half_h)^4;
       1 - ((P_ba{1}(2)-X_b{1}(2))/half_w)^4 - (((P_ba{1}(1)-X_b{1}(1)) + shear*(P_ba{1}(2)-X_b{1}(2)))/half_h)^4;
       1 - ((P_ac{1}(2)-X_a{1}(2))/half_w)^4 - (((P_ac{1}(1)-X_a{1}(1)) + shear*(P_ac{1}(2)-X_a{1}(2)))/half_h)^4;
       1 - ((P_ca{1}(2)-X_c{1}(2))/half_w)^4 - (((P_ca{1}(1)-X_c{1}(1)) + shear*(P_ca{1}(2)-X_c{1}(2)))/half_h)^4;
       1 - ((P_bc{1}(2)-X_b{1}(2))/half_w)^4 - (((P_bc{1}(1)-X_b{1}(1)) + shear*(P_bc{1}(2)-X_b{1}(2)))/half_h)^4;
       1 - ((P_cb{1}(2)-X_c{1}(2))/half_w)^4 - (((P_cb{1}(1)-X_c{1}(1)) + shear*(P_cb{1}(2)-X_c{1}(2)))/half_h)^4];

vars_a = X_a{1};
vars_b = X_b{1};
vars_c = X_c{1};
vars_r = [P_ab{1}; P_ba{1}; 
          P_ac{1}; P_ca{1}; 
          P_bc{1}; P_cb{1}];


p_a = X_a{1};
p_b = X_b{1};
p_c = X_c{1};
p_r = [];

q_a = [];
q_b = [];
q_c = [];


for t = 1:T
    % Introduce symbolic vars.
    U_a{t} = MX.sym(['U_a_' num2str(t)], 2);
    U_b{t} = MX.sym(['U_b_' num2str(t)], 2);
    U_c{t} = MX.sym(['U_c_' num2str(t)], 2);
    X_a{t+1} = MX.sym(['X_a_' num2str(t)], 3);
    X_b{t+1} = MX.sym(['X_b_' num2str(t)], 3);
    X_c{t+1} = MX.sym(['X_c_' num2str(t)], 3);
    
    % Referee primal vars.
    P_ab{t+1} = MX.sym(['P_ab_' num2str(t)], 2);
    P_ba{t+1} = MX.sym(['P_ba_' num2str(t)], 2);
    P_ac{t+1} = MX.sym(['P_ac_' num2str(t)], 2);
    P_ca{t+1} = MX.sym(['P_ca_' num2str(t)], 2);
    P_bc{t+1} = MX.sym(['P_bc_' num2str(t)], 2);
    P_cb{t+1} = MX.sym(['P_cb_' num2str(t)], 2);
   
    vars_a = [vars_a; U_a{t}; X_a{t+1}];
    vars_b = [vars_b; U_b{t}; X_b{t+1}];
    vars_c = [vars_c; U_c{t}; X_c{t+1}];
    vars_r = [vars_r; 
        P_ab{t+1}; P_ba{t+1}; 
        P_ac{t+1}; P_ca{t+1}; 
        P_bc{t+1}; P_cb{t+1}];
    
    % Place vars in lists for easy access later.
    p_a = [p_a; X_a{t+1}];
    p_b = [p_b; X_b{t+1}];
    p_c = [p_c; X_c{t+1}];
    p_r = [p_r; P_ab{t+1}; P_ba{t+1}; P_ac{t+1}; P_ca{t+1}; P_bc{t+1}; P_cb{t+1}];
    q_a = [q_a; U_a{t}];
    q_b = [q_b; U_b{t}];
    q_c = [q_c; U_c{t}];
    
    % Dynamic equations
    X_a_pred = X_a{t} + dt * [X_a{t}(3) + dt/2*U_a{t}(2); ...
                              U_a{t}(1); ...
                              U_a{t}(2)];
    X_b_pred = X_b{t} + dt * [X_b{t}(3) + dt/2*U_b{t}(2); ...
                              U_b{t}(1); ...
                              U_b{t}(2)];
    X_c_pred = X_c{t} + dt * [X_c{t}(3) + dt/2*U_c{t}(2);
                              U_c{t}(1); ...
                              U_c{t}(2)];
    
    % Dynamic constraints.                
    h_a = [h_a; X_a_pred - X_a{t+1}];
    h_b = [h_b; X_b_pred - X_b{t+1}];
    h_c = [h_c; X_c_pred - X_c{t+1}];
    
   
    % Primal variable bounds.
    lb_a = [lb_a; -lat_bound; -a_bound; -inf_const; -inf_const; -inf_const];
    ub_a = [ub_a;  lat_bound; a_bound;   inf_const;  inf_const; inf_const];
    lb_b = [lb_b; -lat_bound; -a_bound; -inf_const; -inf_const; -inf_const];
    ub_b = [ub_b;  lat_bound; a_bound;   inf_const;  inf_const; inf_const];
    lb_c = [lb_c; -lat_bound; -a_bound; -inf_const; -inf_const; -inf_const];
    ub_c = [ub_c;  lat_bound; a_bound;   inf_const;  inf_const; inf_const];
    
    % Collision avoidance constraints.
%     g_a = [g_a; 
%         -2 + ((P_ba{t+1}(2)-X_a{t+1}(2))/half_w)^4 + (((P_ba{t+1}(1)-X_a{t+1}(1)) + shear*(P_ba{t+1}(2)-X_a{t+1}(2)))/half_h)^4;
%         -2 + ((P_ca{t+1}(2)-X_a{t+1}(2))/half_w)^4 + (((P_ca{t+1}(1)-X_a{t+1}(1)) + shear*(P_ca{t+1}(2)-X_a{t+1}(2)))/half_h)^4];
%         -2 + ((P_ba{t+1}(2)-X_a{t+1}(2))/half_w)^4 + (((P_ba{t+1}(1)-X_a{t+1}(1)) + shear*(P_ba{t+1}(2)-X_a{t+1}(2)))/half_h)^4;
%         -2 + ((P_ca{t+1}(2)-X_a{t+1}(2))/half_w)^4 + (((P_ca{t+1}(1)-X_a{t+1}(1)) + shear*(P_ca{t+1}(2)-X_a{t+1}(2)))/half_h)^4];
%         (X_a{t+1}(1:2)-X_b{t+1}(1:2))'*(X_a{t+1}(1:2)-X_b{t+1}(1:2)) - 100;
%         (X_a{t+1}(1:2)-X_c{t+1}(1:2))'*(X_a{t+1}(1:2)-X_c{t+1}(1:2)) - 100];
%         1000*((P_ab{t}-P_ba{t})'*(P_ab{t}-P_ba{t}) - 0.1);
%         1000*((P_ac{t}-P_ca{t})'*(P_ac{t}-P_ca{t}) - 0.1)];
%         1 - ((P_ab{t}(2)-X_a{t+1}(2))/half_w)^4 - (((P_ab{t}(1)-X_a{t+1}(1)) + shear*(P_ab{t}(2)-X_a{t+1}(2)))/half_h)^4;
%         1 - ((P_ac{t}(2)-X_a{t+1}(2))/half_w)^4 - (((P_ac{t}(1)-X_a{t+1}(1)) + shear*(P_ac{t}(2)-X_a{t+1}(2)))/half_h)^4];
%     g_b = [g_b; 
%         -2 + ((P_ab{t+1}(2)-X_b{t+1}(2))/half_w)^4 + (((P_ab{t+1}(1)-X_b{t+1}(1)) + shear*(P_ab{t+1}(2)-X_b{t+1}(2)))/half_h)^4;
%         -2 + ((P_cb{t+1}(2)-X_b{t+1}(2))/half_w)^4 + (((P_cb{t+1}(1)-X_b{t+1}(1)) + shear*(P_cb{t+1}(2)-X_b{t+1}(2)))/half_h)^4];
%         -2 + ((P_ab{t+1}(2)-X_b{t+1}(2))/half_w)^4 + (((P_ab{t+1}(1)-X_b{t+1}(1)) + shear*(P_ab{t+1}(2)-X_b{t+1}(2)))/half_h)^4;
%         -2 + ((P_cb{t+1}(2)-X_b{t+1}(2))/half_w)^4 + (((P_cb{t+1}(1)-X_b{t+1}(1)) + shear*(P_cb{t+1}(2)-X_b{t+1}(2)))/half_h)^4];
%         (X_b{t+1}(1:2)-X_a{t+1}(1:2))'*(X_b{t+1}(1:2)-X_a{t+1}(1:2)) - 100;
%         (X_b{t+1}(1:2)-X_c{t+1}(1:2))'*(X_b{t+1}(1:2)-X_c{t+1}(1:2)) - 100];
%         1*((X_b{t+1}(1:2)-P_ab{t})'*(X_b{t+1}(1:2)-P_ab{t}) - (X_b{t+1}(1:2)-P_ba{t})'*(X_b{t+1}(1:2)-P_ba{t}) - 0.25);
%         1*((X_b{t+1}(1:2)-P_cb{t})'*(X_b{t+1}(1:2)-P_cb{t}) - (X_b{t+1}(1:2)-P_bc{t})'*(X_b{t+1}(1:2)-P_bc{t}) - 0.25)];
%         1000*((P_ba{t}-P_ab{t})'*(P_ba{t}-P_ab{t}) - 0.1);
%         1000*((P_bc{t}-P_cb{t})'*(P_bc{t}-P_cb{t}) - 0.1)];
%         1 - ((P_ba{t}(2)-X_b{t+1}(2))/half_w)^4 - (((P_ba{t}(1)-X_b{t+1}(1)) + shear*(P_ba{t}(2)-X_b{t+1}(2)))/half_h)^4;
%         1 - ((P_bc{t}(2)-X_b{t+1}(2))/half_w)^4 - (((P_bc{t}(1)-X_b{t+1}(1)) + shear*(P_bc{t}(2)-X_b{t+1}(2)))/half_h)^4];
%     g_c = [g_c;
%         1 - ((P_ca{t}(2)-X_c{t+1}(2))/half_w)^4 - (((P_ca{t}(1)-X_c{t+1}(1)) + shear*(P_ca{t}(2)-X_c{t+1}(2)))/half_h)^4;
%         1 - ((P_cb{t}(2)-X_c{t+1}(2))/half_w)^4 - (((P_cb{t}(1)-X_c{t+1}(1)) + shear*(P_cb{t}(2)-X_c{t+1}(2)))/half_h)^4];
    % Referee constraints. (Solve closest point problem).
    g_r = [g_r; 
        1 - ((P_ab{t+1}(2)-X_a{t+1}(2))/half_w)^4 - (((P_ab{t+1}(1)-X_a{t+1}(1)) + shear*(P_ab{t+1}(2)-X_a{t+1}(2)))/half_h)^4;
        1 - ((P_ba{t+1}(2)-X_b{t+1}(2))/half_w)^4 - (((P_ba{t+1}(1)-X_b{t+1}(1)) + shear*(P_ba{t+1}(2)-X_b{t+1}(2)))/half_h)^4;
        1 - ((P_ac{t+1}(2)-X_a{t+1}(2))/half_w)^4 - (((P_ac{t+1}(1)-X_a{t+1}(1)) + shear*(P_ac{t+1}(2)-X_a{t+1}(2)))/half_h)^4;
        1 - ((P_ca{t+1}(2)-X_c{t+1}(2))/half_w)^4 - (((P_ca{t+1}(1)-X_c{t+1}(1)) + shear*(P_ca{t+1}(2)-X_c{t+1}(2)))/half_h)^4;
        1 - ((P_bc{t+1}(2)-X_b{t+1}(2))/half_w)^4 - (((P_bc{t+1}(1)-X_b{t+1}(1)) + shear*(P_bc{t+1}(2)-X_b{t+1}(2)))/half_h)^4;
        1 - ((P_cb{t+1}(2)-X_c{t+1}(2))/half_w)^4 - (((P_cb{t+1}(1)-X_c{t+1}(1)) + shear*(P_cb{t+1}(2)-X_c{t+1}(2)))/half_h)^4];
    
    % Cost function design.
    f_a = f_a + ...
        0.01*(U_a{t}(1)'*U_a{t}(1)) + ... % Minimize lateral velocity
        0.01*(U_a{t}(2)'*U_a{t}(2)) + ... % Minimize acceleration
        0.01*(X_a{t+1}(2) - 0)^2 * (X_a{t+1}(2)-5)^2 + ... % Stay in lane centers
        .01*(X_a{t+1}(3) - va)^2 + ... % Go the desired speed 
        2*(X_b{t+1}(3) - vb)^2; % Let pb go the desired speed
    f_b = f_b + ...
        0.01*(U_b{t}(1)'*U_b{t}(1)) + ...
        0.01*(U_b{t}(2)'*U_b{t}(2)) + ...
        0.01*(X_b{t+1}(2) - 0)^2*(X_b{t+1}(2)-5)^2 + ... % Stay in lane centers
        .01*(X_b{t+1}(3) - vb)^2 + ... 
        2*(X_a{t+1}(3) - va)^2;
    f_c = f_c + ...
        0.1*(U_c{t}'*U_c{t}) + ...
        0.1*(-X_c{t+1}(1)) + ... % Forward progress
        2*(X_c{t+1}(3) -vc)^2 + ...
        0.01*(X_c{t+1}(2) - 0)^2*(X_c{t+1}(2)-5)^2;
    
    f_r = f_r + ...
        1*(P_ab{t+1}-P_ba{t+1})'*(P_ab{t+1}-P_ba{t+1}) + ...
        1*(P_ac{t+1}-P_ca{t+1})'*(P_ac{t+1}-P_ca{t+1}) + ...
        1*(P_bc{t+1}-P_cb{t+1})'*(P_bc{t+1}-P_cb{t+1});
end

lam_a = MX.sym('lam_a', size(h_a,1));
lam_b = MX.sym('lam_b', size(h_b,1));
lam_c = MX.sym('lam_c', size(h_c,1));
lam_r = MX.sym('lam_r', size(h_r,1));

mu_a = MX.sym('mu_a', size(g_a,1));
mu_b = MX.sym('mu_b', size(g_b,1));
mu_c = MX.sym('mu_c', size(g_c,1));
mu_r = MX.sym('mu_r', size(g_r,1));

Lag_a = f_a; 
Lag_b = f_b;
Lag_c = f_c;
Lag_r = f_r;

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

if size(h_c,1) > 0
    Lag_c = Lag_c - lam_c'*h_c;
end
if size(g_c,1) > 0
    Lag_c = Lag_c - mu_c'*g_c;
end

if size(h_r,1) > 0
    Lag_r = Lag_r - lam_r'*h_r;
end
if size(g_r,1) > 0
    Lag_r = Lag_r - mu_r'*g_r;
end



gLa = gradient(Lag_a, vars_a);
gLb = gradient(Lag_b, vars_b);
gLc = gradient(Lag_c, vars_c);
gLr = gradient(Lag_r, vars_r);

gL = [gLa; gLb; gLc; gLr; h_a; h_b; h_c; h_r; g_a; g_b; g_c; g_r];

all_vars = [vars_a; vars_b; vars_c; vars_r; lam_a; lam_b; lam_c; lam_r; mu_a; mu_b; mu_c; mu_r];
HH = jacobian(gL, all_vars);
% spy(HH)

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});
evalFa = Function('obj_a', {all_vars}, {f_a});
evalFb = Function('obj_b', {all_vars}, {f_b});
evalha = Function('const_a', {all_vars}, {h_a});
evalhb = Function('const_b', {all_vars}, {h_b});
evalhc = Function('const_c', {all_vars}, {h_c});
evalhr = Function('const_r', {all_vars}, {h_r});

evalga = Function('const_a', {all_vars}, {g_a});
evalgb = Function('const_b', {all_vars}, {g_b});
evalgc = Function('const_c', {all_vars}, {g_c});
evalgr = Function('const_r', {all_vars}, {g_r});

extract_pa = Function('extract_pa', {all_vars}, {p_a});
extract_pb = Function('extract_pb', {all_vars}, {p_b});
extract_pc = Function('extract_pc', {all_vars}, {p_c});
extract_pr = Function('extract_pr', {all_vars}, {vars_r});

extract_qa = Function('extract_qa', {all_vars}, {q_a});
extract_qb = Function('extract_qb', {all_vars}, {q_b});
extract_qc = Function('extract_qc', {all_vars}, {q_c});

extract_lam_a = Function('extract_mu_a', {all_vars}, {lam_a});
extract_lam_b = Function('extract_mu_b', {all_vars}, {lam_b});
extract_lam_c = Function('extract_mu_c', {all_vars}, {lam_c});

extract_mu_a = Function('extract_mu_a', {all_vars}, {mu_a});
extract_mu_b = Function('extract_mu_b', {all_vars}, {mu_b});
extract_mu_c = Function('extract_mu_c', {all_vars}, {mu_c});


ones_a = ones(size(vars_a,1),1);
ones_lam_a = ones(size(lam_a,1),1);
ones_mu_a = ones(size(mu_a,1),1);

ones_b = ones(size(vars_b,1),1);
ones_lam_b = ones(size(lam_b,1),1);
ones_mu_b = ones(size(mu_b,1),1);

ones_c = ones(size(vars_c,1),1);
ones_lam_c = ones(size(lam_c,1),1);
ones_mu_c = ones(size(mu_c,1),1);

ones_r = ones(size(vars_r,1),1);
ones_lam_r = ones(size(lam_r,1),1);
ones_mu_r = ones(size(mu_r,1),1);

lo = [lb_a; lb_b; lb_c; -inf_const*ones_r;
    -inf_const*ones_lam_a;
    -inf_const*ones_lam_b;
    -inf_const*ones_lam_c;
    -inf_const*ones_lam_r;
    0*ones_mu_a;
    0*ones_mu_b;
    0*ones_mu_c;
    0*ones_mu_r];
hi = [ub_a; ub_b; ub_c; inf_const*ones_r;
    inf_const*ones_lam_a;
    inf_const*ones_lam_b;
    inf_const*ones_lam_c;
    inf_const*ones_lam_r;
    inf_const*ones_mu_a;
    inf_const*ones_mu_b;
    inf_const*ones_mu_c;
    inf_const*ones_mu_r];

z0 = zeros(size(all_vars,1),1);
xx_a = x0_a;
xx_b = x0_b;
xx_c = x0_c;
init_a = x0_a;
init_b = x0_b;
init_c = x0_c;
init_r = [xx_a(1:2);xx_b(1:2);xx_a(1:2);xx_c(1:2);xx_b(1:2);xx_c(1:2)];

for t = 1:T
    % Stay in lane params.
    u_a = [0;0];
   
    if (t < 16)
        u_b = [0;0];
    else
        u_b = [0;0];
    end
    u_c = [0;0];
    xx_a = xx_a + dt * [xx_a(3)+dt/2*u_a(2);u_a];
    init_a = [init_a; u_a; xx_a];
    xx_b = xx_b + dt * [xx_b(3)+dt/2*u_b(2);u_b];
    init_b = [init_b; u_b; xx_b];
    xx_c = xx_c + dt * [xx_c(3)+dt/2*u_c(2);u_c];
    init_c = [init_c; u_c; xx_c];
    init_r = [init_r; xx_a(1:2);xx_b(1:2);xx_a(1:2);xx_c(1:2);xx_b(1:2);xx_c(1:2)];
end

XA = init_a(1:5:end);
YA = init_a(2:5:end);

XB = init_b(1:5:end);
YB = init_b(2:5:end);

XC = init_c(1:5:end);
YC = init_c(2:5:end);

S = size(vars_a,1);
z0(1:S) = init_a;
z0(S+1:2*S) = init_b;
z0(2*S+1:3*S) = init_c;
z0(3*S+1:3*S+size(vars_r,1)) = init_r;

n = size(z0,1);
m = size([vars_a; vars_b; vars_c; vars_r; lam_a; lam_b; lam_c; lam_r],1);
nnz = nnz(evalH(z0));

init_loc_a = 1;
init_loc_b = 1;
dim_state = 0;
dim_control = 0;
num_iterations=1;
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
ZZ = readz('/Users/forrest/code/packages/pathlib-master/examples/C/z.csv');

for iter = 1:1
    z = ZZ((iter-1)*n+1:iter*n);
    

%     ineq_violation_a = min(0,min(evalga(z)))
%     ineq_violation_b = min(0,min(evalgb(z)))
%     ineq_violation_c = min(0,min(evalgc(z)))
%     ineq_violation_r = min(0,min(evalgr(z)))
% 
%     eq_violation_a = max(abs(evalha(z)))
%     eq_violation_b = max(abs(evalhb(z)))
%     eq_violation_c = max(abs(evalhc(z)))
%     eq_violation_r = max(abs(evalhr(z)))

    PA = full(extract_pa(z));
    PB = full(extract_pb(z));
    PC = full(extract_pc(z));
    PR = full(extract_pr(z));

    QA = full(extract_qa(z));
    QB = full(extract_qb(z));
    QC = full(extract_qc(z));

    XA = PA(1:3:end);
    YA = PA(2:3:end);
    VA = PA(3:3:end);
    AA = QA(2:2:end);
    LA = QA(1:2:end);

    XB = PB(1:3:end);
    YB = PB(2:3:end);
    VB = PB(3:3:end);
    AB = QB(2:2:end);
    LB = QB(1:2:end);

    XC = PC(1:3:end);
    YC = PC(2:3:end);
    VC = PC(3:3:end);
    AC = QC(2:2:end);
    LC = QC(1:2:end);

    P_AB_X = PR(1:12:end); P_AB_Y = PR(2:12:end);
    P_BA_X = PR(3:12:end); P_BA_Y = PR(4:12:end);
    P_AC_X = PR(5:12:end); P_AC_Y = PR(6:12:end);
    P_CA_X = PR(7:12:end); P_CA_Y = PR(8:12:end);
    P_BC_X = PR(9:12:end); P_BC_Y = PR(10:12:end);
    P_CB_X = PR(11:12:end); P_CB_Y = PR(12:12:end);

    figure; hold on;
    plot([-2.5,-2.5],[-5,500],'k');
    plot([2.5,2.5],[-5,500],'k');
    plot([7.5,7.5],[-5,500],'k');
    for t = 1:T+1
        hold on;
        path_a = plot(YA(t,1),XA(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        int_ab = plot(P_AB_Y(t),P_AB_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        int_ac = plot(P_AC_Y(t),P_AC_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
        [Xsrc_a, Ysrc_a] = drawShearedRecircle(XA(t,1),YA(t,1),shear,half_w,half_h);
        ca = plot(Xsrc_a,Ysrc_a,'b');
    %     ca1 = viscircles([YA(t,1),XA(t,1)]+off_for,bounding_radius,'Color', 'b');
    %     ca2 = viscircles([YA(t,1),XA(t,1)]+off_bac,bounding_radius,'Color', 'b');
        path_b = plot(YB(t,1),XB(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        int_ba = plot(P_BA_Y(t),P_BA_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        int_bc = plot(P_BC_Y(t),P_BC_X(t), 'or','MarkerSize', 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
        [Xsrc_b, Ysrc_b] = drawShearedRecircle(XB(t,1),YB(t,1),shear,half_w,half_h);
        cb = plot(Xsrc_b,Ysrc_b,'r');
    %     cb1 = viscircles([YB(t,1),XB(t,1)]+off_for,bounding_radius,'Color', 'r');
    %     cb2 = viscircles([YB(t,1),XB(t,1)]+off_bac,bounding_radius,'Color', 'r');
        path_c = plot(YC(t,1),XC(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
        int_ca = plot(P_CA_Y(t),P_CA_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        int_cb = plot(P_CB_Y(t),P_CB_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        [Xsrc_c, Ysrc_c] = drawShearedRecircle(XC(t,1),YC(t,1),shear,half_w,half_h);
        cc = plot(Xsrc_c,Ysrc_c,'g');
    %     cc1 = viscircles([YC(t,1),XC(t,1)]+off_for,bounding_radius,'Color', 'g');
    %     cc2 = viscircles([YC(t,1),XC(t,1)]+off_bac,bounding_radius,'Color', 'g');
        pos = XB(t,1);
        axis([-30 30 pos-30 pos+30]);
        pause(dt);
        delete(ca);
        delete(cb);
        delete(cc);
        delete(path_a);
        delete(path_b);
        delete(path_c);
        delete(int_ab);
        delete(int_ac);
        delete(int_ba);
        delete(int_bc);
        delete(int_ca);
        delete(int_cb);
    %     delete(ca2);
    %     delete(cb2);
    %     delete(cc2);
    end

end

