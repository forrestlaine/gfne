clear; clc; close all;
import casadi.*

inf_const = 999999999999.9;
closed_loop = false;

dt = 0.25;
T = 1;
num_iterations = 1;
S = 40;

x0_a = [5;0;-1;0];
x0_b = [0;5;0;0];

% f_a = 0;
% f_b = 0;
% 
% h_a = [];
% h_b = [];

p_a = x0_a;
p_b = x0_b;

q_a = [];
q_b = [];

for s = 1:S
    vars_a{s} = [];
    vars_b{s} = [];
    vars_only_a{s} = [];
    vars_only_b{s} = [];
    
    f_a{s} = 0;
    f_b{s} = 0;
    h_a{s} = [];
    h_b{s} = [];
    h_c{s} = [];
    g_a{s} = [];
    g_b{s} = [];
    
    for t = 1:T
        
        U_a{s,t} = MX.sym(['U_a_' num2str(s) '_' num2str(t)], 2);
        U_b{s,t} = MX.sym(['U_b_' num2str(s) '_' num2str(t)], 2);

        X_a{s,t} = MX.sym(['X_a_' num2str(s) '_' num2str(t)], 4);
        X_b{s,t} = MX.sym(['X_b_' num2str(s) '_' num2str(t)], 4);

        vars_a{s} = [vars_a{s}; U_a{s,t}; X_a{s,t}; X_b{s,t}];
        vars_b{s} = [vars_b{s}; U_b{s,t}; X_a{s,t}; X_b{s,t}];
        
        vars_only_a{s} = [vars_only_a{s}; U_a{s,t}; X_a{s,t}];
        vars_only_b{s} = [vars_only_b{s}; U_b{s,t}; X_b{s,t}];


        % Place vars in lists for easy access later.
        p_a = [p_a; X_a{s,t}];
        p_b = [p_b; X_b{s,t}];

        q_a = [q_a; U_a{s,t}];
        q_b = [q_b; U_b{s,t}];

        if (t==1)
            if (s==1)
                X_a_prev = x0_a;
                X_b_prev = x0_b;
            else
                X_a_prev = X_a{s-1,T};
                X_b_prev = X_b{s-1,T};
            end
        else
            X_a_prev = X_a{s,t-1};
            X_b_prev = X_b{s,t-1};
        end
        
        dubins_dyn_a = [X_a_prev(3)*cos(X_a_prev(4)); 
                        X_a_prev(3)*sin(X_a_prev(4));
                        U_a{s,t}(1); 
                        U_a{s,t}(2)];
        dubins_dyn_b = [X_b_prev(3)*cos(X_b_prev(4)); 
                        X_b_prev(3)*sin(X_b_prev(4)); 
                        U_b{s,t}(1);
                        U_b{s,t}(2)];
        pointmass_dyn_a = [X_a_prev(3); 
                        X_a_prev(4);
                        U_a{s,t}(1); 
                        U_a{s,t}(2)];
        pointmass_dyn_b = [X_b_prev(3); 
                        X_b_prev(4); 
                        U_b{s,t}(1);
                        U_b{s,t}(2)];
                    
%         X_a_pred = X_a_prev + dt * dubins_dyn_a;
%         X_b_pred = X_b_prev + dt * dubins_dyn_b;
        X_a_pred = X_a_prev + dt * pointmass_dyn_a;
        X_b_pred = X_b_prev + dt * pointmass_dyn_b;

        % Dynamic constraints.                
        h_a{s} = [h_a{s}; X_a_pred - X_a{s,t}];
        h_b{s} = [h_b{s}; X_b_pred - X_b{s,t}];
        h_c{s} = [h_c{s}; X_a_pred - X_a{s,t}; X_b_pred - X_b{s,t}];
        
        
        % Cost function design.
        f_a{s} = f_a{s} + ...
            1*(U_a{s,t}'*U_a{s,t}) + ...
            1*((X_b{s,t}(1:2))'*(X_b{s,t}(1:2)));
        f_b{s} = f_b{s} + ...
            1*(U_b{s,t}'*U_b{s,t}) + ...
             1*((X_b{s,t}(1:2) - X_a{s,t}(1:2))' * (X_b{s,t}(1:2) - X_a{s,t}(1:2)));

    end
end

h_a_all = [];
h_b_all = [];
h_c_all = [];
for s = 1:S
    h_a_all = [h_a_all; h_a{s}];
    h_b_all = [h_b_all; h_b{s}];
    h_c_all = [h_c_all; h_c{s}];
end

lam_a_all = MX.sym('lam_a', size(h_c_all,1));
lam_b_all = MX.sym('lam_a', size(h_c_all,1));

for s = 1:S
    if (s < S)
        g_a{s} = h_a{s+1}(1:4);
        g_b{s} = h_b{s+1}(1:4);
        g_c{s} = h_c{s+1}(1:8);
        gam_a{s} = MX.sym(['gam_a_' num2str(s)], 8);
        gam_b{s} = MX.sym(['gam_b_' num2str(s)], 8);
    end
    lam_a{s} = MX.sym(['lam_a_' num2str(s)], size(h_c{s},1));
    lam_b{s} = MX.sym(['lam_b_' num2str(s)], size(h_c{s},1));
end

% Set up necessary condition vector 
for s = 1:S
    Lag_a{s} = f_a{s} - lam_a{s}'*h_c{s};
    Lag_b{s} = f_b{s} - lam_b{s}'*h_c{s};
    
    if (s < S)
        Lag_a{s} = Lag_a{s} - gam_a{s}'*g_c{s};
        Lag_b{s} = Lag_b{s} - gam_b{s}'*g_c{s};
        infl_a{s} = gradient(gam_a{s}'*g_c{s}, vars_a{s});
        infl_b{s} = gradient(gam_b{s}'*g_c{s}, vars_b{s});
    end
    
    gLa{s} = gradient(Lag_a{s},vars_a{s});
    gLb{s} = gradient(Lag_b{s},vars_b{s});
end

gL = [];

for s = 1:S
    pL{s} = [gradient(f_a{s} - lam_a{s}'*h_c{s}, vars_a{s});
             gradient(f_b{s} - lam_b{s}'*h_c{s}, vars_b{s});
             h_c{s}];
    sL{s} = [gradient(f_b{s} - lam_b{s}'*h_c{s}, vars_a{s});
             gradient(f_a{s} - lam_a{s}'*h_c{s}, vars_b{s});
             h_a{s};
             h_b{s}];
    stage_vars{s} = [vars_only_a{s}; vars_only_b{s}; lam_a{s}; lam_b{s}];
    sys{s} = jacobian(pL{s},stage_vars{s});
    sys_inv{s} = inv(sys{s});
%     s_sys_inv{s} = inv(jacobian(sL{s},stage_vars{s}));
    na = size(vars_only_a{s},1);
    nb = size(vars_only_b{s},1);
    
    ma = size(vars_a{s},1);
    mb = size(vars_b{s},1);
    m = 2;
    l = 8;
%     na = size(vars_only_a{s},1);
%     nb = size(vars_b{s},1);
%     la = size(lam_a{s},1);
%     lb = size(lam_b{s},1);
    d_vars_ua_dh{s} = sys_inv{s}(1:2, ma+mb+1:end);
    d_vars_xa_dh{s} = sys_inv{s}(3:6, ma+mb+1:end);
    d_vars_ub_dh{s} = sys_inv{s}(7:8, ma+mb+1:end);
    d_vars_xb_dh{s} = sys_inv{s}(9:12, ma+mb+1:end);
    
    d_vars_aa_dh{s} = [d_vars_ua_dh{s}; d_vars_xa_dh{s}; d_vars_xb_dh{s}];
    d_vars_bb_dh{s} = [d_vars_ub_dh{s}; d_vars_xa_dh{s}; d_vars_xb_dh{s}];
     
    d_vars_all_dh{s} = [d_vars_ua_dh{s}; d_vars_xa_dh{s}; d_vars_ub_dh{s}; d_vars_xb_dh{s}];
    
%     d_vars_aa_dh{s} = sys_inv{s}(1:na, ma+mb+1:end);
%     d_vars_bb_dh{s} = sys_inv{s}(na+1:na+nb, ma+mb+1:end);
%     d_vars_ab_dh{s} = sys_inv{s}(1:na, na+nb+la+1:na+nb+la+lb);
%     d_vars_ba_dh{s} = sys_inv{s}(na+1:na+nb, na+nb+1:na+nb+la);
%     d_vars_bb_dh{s} = sys_inv{s}(na+1:na+nb, na+nb+la+1:na+nb+la+lb);
end

for s = 1:S
    gL = [gL; gLa{s}];
end
for s = 1:S
    gL = [gL; gLb{s}];
end
for s = 1:S
    gL = [gL; h_c{s}];
end
% for s = 1:S
%     gL = [gL; h_b{s}];
% end
% first_stage_conditions = gL;
% for s = 1:S-1
%     first_stage_conditions = [first_stage_conditions; gam_a{s}; gam_b{s}];
% end


% Set up variable vector (matching condition order)
all_vars = [];
for s = 1:S
    all_vars = [all_vars; vars_only_a{s}];
end
for s = 1:S
    all_vars = [all_vars; vars_only_b{s}];
end
traj_size = size(all_vars,1);
for s = 1:S
    all_vars = [all_vars; lam_a{s}];
end
for s = 1:S
    all_vars = [all_vars; lam_b{s}];
end
traj_lam_size = size(all_vars,1);
for s = 1:S-1
    all_vars = [all_vars; gam_a{s}; gam_b{s}];
end

% preHH = jacobian(first_stage_conditions, all_vars);
% preHHi = inv(preHH);
ind = 0;
jnd = traj_size;
for s = 1:S
    n = size(vars_a{s},1);
%     d_vars_a_dh{s} = preHHi(ind+1:ind+n,jnd+1:jnd+4);
    if (s < S)
        func = f_a{s}-gam_a{s}'*g_c{s};
    else
        func = f_a{s};
    end
    traj_grad_a{s} = gradient(func,[vars_only_a{s}; vars_only_b{s}]);
    shadow_grad_b{s} = gradient(func,vars_b{s});

%     ind = ind+n;
%     jnd = jnd+size(h_a{s},1);
end
for s = 1:S
    n = size(vars_b{s},1);
%     d_vars_b_dh{s} = preHHi(ind+1:ind+n,jnd+1:jnd+4);
    if (s < S)
        func = f_b{s}-gam_b{s}'*g_c{s};
    else
        func = f_b{s};
    end
    traj_grad_b{s} = gradient(func,[vars_only_a{s}; vars_only_b{s}]);
    shadow_grad_a{s} = gradient(func,vars_a{s}); % change in cost_b wrt traj_a
%     ind = ind+n;
%     jnd = jnd+size(h_b{s},1);
end

for s = 1:S-1
    if closed_loop
        gL = [gL; 
            gam_a{s}-d_vars_aa_dh{s+1}'*traj_grad_a{s+1}-d_vars_ba_dh{s+1}'*shadow_grad_b{s+1}; 
            gam_b{s}-d_vars_bb_dh{s+1}'*traj_grad_b{s+1}-d_vars_ab_dh{s+1}'*shadow_grad_a{s+1}];
    else
        gL = [gL; 
            gam_a{s}-d_vars_all_dh{s+1}'*traj_grad_a{s+1}; 
            gam_b{s}-d_vars_all_dh{s+1}'*traj_grad_b{s+1}];        
    end
end

HH = jacobian(gL, all_vars);
figure;
spy(HH);

% evaldi = Function('d1', {all_vars}, {preHHi});
% z = ones(size(all_vars));
evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});

extract_pa = Function('extract_pa', {all_vars}, {p_a});
extract_pb = Function('extract_pb', {all_vars}, {p_b});

extract_qa = Function('extract_qa', {all_vars}, {q_a});
extract_qb = Function('extract_qb', {all_vars}, {q_b});
for s = 1:S
%     get_pL{s} = Function(['extract_pL_' num2str(s)], {all_vars}, {sys_inv{s}});
    dvars_all{s} = Function(['extract_dvars_ab_' num2str(s)], {all_vars}, {d_vars_all_dh{s}});
    if s < S
        get_gam_a{s} = Function('get_gam_a', {all_vars}, {gam_a{s}});
        get_gam_b{s} = Function('get_gam_b', {all_vars}, {gam_b{s}});
        get_g{s} = Function('get_g', {all_vars}, {g_c{s}});
        get_infl_a{s} = Function('get_infl_a', {all_vars}, {infl_a{s}});
        get_infl_b{s} = Function('get_infl_b', {all_vars}, {infl_b{s}});
        get_sys{s} = Function('get_sys', {all_vars}, {sys_inv{s}});
    end
%     dvars_b{s} = Function(['extract_dvars_ba_' num2str(s)], {all_vars}, {d_vars_bb_dh{s}});
    tgrad_a{s} = Function(['extract_shadow_a_' num2str(s)], {all_vars}, {traj_grad_a{s}});
    tgrad_b{s} = Function(['extract_shadow_b_' num2str(s)], {all_vars}, {traj_grad_b{s}});
end


extract_qb = Function('extract_qb', {all_vars}, {q_b});

lo = -inf_const*ones(size(all_vars));
hi = inf_const*ones(size(all_vars));
z0 = zeros(size(all_vars));

% Specs needed by C program
n = size(z0,1);
m = n;
nnz = nnz(evalH(z0));
init_loc_a = 0;
init_loc_b = size(vars_a,1);
dim_state = size(X_a{1},1);
dim_control = size(U_a{1},1);
specs = [n;m;nnz;init_loc_a;init_loc_b;dim_state;dim_control;T;num_iterations];

% Write to files so Path solver can operate
writez(lo, '/Users/forrest/code/packages/pathlib-master/examples/C/lo.txt');
writez(hi, '/Users/forrest/code/packages/pathlib-master/examples/C/hi.txt');
writez(z0, '/Users/forrest/code/packages/pathlib-master/examples/C/z0.txt');
write_nums(specs, '/Users/forrest/code/packages/pathlib-master/examples/C/spec.txt');

% Compile C code
opts = struct('with_header', true, 'casadi_int', 'int');
C = CodeGenerator('gen.c', opts);
C.add(evalG);
C.add(evalH);
C.generate();

disp 'Generated game!';

% Run path solver
cmd = "cd /Users/forrest/code/packages/pathlib-master/examples/C && " + ...
      "export DYLD_LIBRARY_PATH=/Users/forrest/code/packages/pathlib-master/lib/osx && " + ...
      "export PATH_LICENSE_STRING='2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0' && " + ...
      "cp ~/Documents/MATLAB/fb_games/test_single_player/gen.c . && make && ./game_shared";
[status,cmdout] = system(cmd)
%%
ZZ = readz('/Users/forrest/code/packages/pathlib-master/examples/C/z.csv');

figure; hold on;
% viscircles([0,0],2);
for iter = 1:num_iterations
    z = ZZ((iter-1)*n+1:iter*n);
%     ineq_violation_a = min(0,min(evalga(z)))
%     ineq_violation_b = min(0,min(evalgb(z)))
% 
%     eq_violation_a = max(abs(evalha(z)))
%     eq_violation_b = max(abs(evalhb(z)))

    PA = full(extract_pa(z));
    PB = full(extract_pb(z));

    QA = full(extract_qa(z));
    QB = full(extract_qb(z));
%     extract_lam_a(z) - extract_shlam_a(z)
%     extract_lam_b(z) - extract_shlam_b(z)

    XA = PA(1:4:end);
    YA = PA(2:4:end);
    AA = QA(2:2:end);
    LA = QA(1:2:end);

    XB = PB(1:4:end);
    YB = PB(2:4:end);
    AB = QB(2:2:end);
    LB = QB(1:2:end);
    
    [XA YA XB YB]

    % P_AB_X = PR(1:12:end); P_AB_Y = PR(2:12:end);
    % P_BA_X = PR(3:12:end); P_BA_Y = PR(4:12:end);
    % P_AC_X = PR(5:12:end); P_AC_Y = PR(6:12:end);
    % P_CA_X = PR(7:12:end); P_CA_Y = PR(8:12:end);
    % P_BC_X = PR(9:12:end); P_BC_Y = PR(10:12:end);
    % P_CB_X = PR(11:12:end); P_CB_Y = PR(12:12:end);

    % plot([-2.5,-2.5],[-5,500],'k');
    % plot([2.5,2.5],[-5,500],'k');
    % plot([7.5,7.5],[-5,500],'k');
   
    path_a = plot(XA(:,1), YA(:,1),'or', 'MarkerSize', 1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
    path_b = plot(XB(:,1),YB(:,1),'or', 'MarkerSize', 1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'r');

%     spot_a = plot(XA(1,1), YA(1,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     spot_b = plot(XB(1,1),YB(1,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     cur_a = viscircles([XA(1,1),YA(1,1)],1,'Color', 'b');
%     cur_b = viscircles([XB(1,1),YB(1,1)],1,'Color', 'r');
    time = text(XA(1,1)+15, YA(1,1)+15,string(iter*dt));
    axis([XA(1,1)-20 XA(1,1)+20 YA(1,1)-20 YA(1,1)+20]);
%     pause(dt);
%     delete(path_a);
%     delete(path_b);
%     delete(cur_a);
%     delete(cur_b);
%     delete(time);
%     for t = 1:T+1
%         hold on;
%         path_a = plot(XA(t,1),YA(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     %     int_ab = plot(P_AB_Y(t),P_AB_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     %     int_ac = plot(P_AC_Y(t),P_AC_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%     %     [Xsrc_a, Ysrc_a] = drawShearedRecircle(XA(t,1),YA(t,1),shear,half_w,half_h);
%     %     ca = plot(Xsrc_a,Ysrc_a,'b');
%     %     ca1 = viscircles([YA(t,1),XA(t,1)]+off_for,bounding_radius,'Color', 'b');
%     %     ca2 = viscircles([YA(t,1),XA(t,1)]+off_bac,bounding_radius,'Color', 'b');
%         path_b = plot(XB(t,1),YB(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     %     int_ba = plot(P_BA_Y(t),P_BA_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     %     int_bc = plot(P_BC_Y(t),P_BC_X(t), 'or','MarkerSize', 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%     %     [Xsrc_b, Ysrc_b] = drawShearedRecircle(XB(t,1),YB(t,1),shear,half_w,half_h);
%     %     cb = plot(Xsrc_b,Ysrc_b,'r');
%     %     cb1 = viscircles([YB(t,1),XB(t,1)]+off_for,bounding_radius,'Color', 'r');
%     %     cb2 = viscircles([YB(t,1),XB(t,1)]+off_bac,bounding_radius,'Color', 'r');
%     %     path_c = plot(YC(t,1),XC(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%     %     int_ca = plot(P_CA_Y(t),P_CA_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     %     int_cb = plot(P_CB_Y(t),P_CB_X(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     %     [Xsrc_c, Ysrc_c] = drawShearedRecircle(XC(t,1),YC(t,1),shear,half_w,half_h);
%     %     cc = plot(Xsrc_c,Ysrc_c,'g');
%     %     cc1 = viscircles([YC(t,1),XC(t,1)]+off_for,bounding_radius,'Color', 'g');
%     %     cc2 = viscircles([YC(t,1),XC(t,1)]+off_bac,bounding_radius,'Color', 'g');
%     %     pos = XB(t,1);
%         axis([-30 30 -30 30]);
%         pause(dt);
%     %     delete(ca);
%     %     delete(cb);
%     %     delete(cc);
%         delete(path_a);
%         delete(path_b);
%     %     delete(path_c);
%     %     delete(int_ab);
%     %     delete(int_ac);
%     %     delete(int_ba);
%     %     delete(int_bc);
%     %     delete(int_ca);
%     %     delete(int_cb);
%     %     delete(ca2);
%     %     delete(cb2);
%     %     delete(cc2);
%     end
end


