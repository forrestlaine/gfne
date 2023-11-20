clear; clc; close all;
import casadi.*

inf_const = 999999999999.9;

dt = 0.25;
T = 40;
num_iterations = 1;

X_a{1} = MX.sym(['X_a_' num2str(0)], 4);
X_b{1} = MX.sym(['X_b_' num2str(0)], 4);


% Referee primal vars.
% P_ab{1} = MX.sym(['P_ab_' num2str(0)], 2);
% P_ba{1} = MX.sym(['P_ba_' num2str(0)], 2);
% P_ac{1} = MX.sym(['P_ac_' num2str(0)], 2);
% P_ca{1} = MX.sym(['P_ca_' num2str(0)], 2);
% P_bc{1} = MX.sym(['P_bc_' num2str(0)], 2);
% P_cb{1} = MX.sym(['P_cb_' num2str(0)], 2);

% 
% half_w = 1.5;
% half_h = 5.0;
% shear = 0.2;

% vehicle_width = 2.75;
% vehicle_length = 3.5;
% 
% offset_forward = [vehicle_length/2; -vehicle_width/8];
% offset_backward = [-vehicle_length/2; vehicle_width/8];
% bounding_radius = sqrt((5/8*vehicle_width)^2);
% min_dist_sq = (2*bounding_radius)^2;

% min_dist_sq = 20;
% lat_bound = inf_const; % m/s
a_bound_a = 1; % m/s^2
a_bound_b = 1;
v_bound = inf_const;

x0_a = [5;0;-1;0];
x0_b = [0;5;0;0];

f_a = 0;
f_b = 0;

% f_r = 1*(P_ab{1}-P_ba{1})'*(P_ab{1}-P_ba{1}) + ...
%       1*(P_ac{1}-P_ca{1})'*(P_ac{1}-P_ca{1}) + ...
%       1*(P_bc{1}-P_cb{1})'*(P_bc{1}-P_cb{1});

h_a = X_a{1}-x0_a;
h_b = X_b{1}-x0_b;
% h_c = [];
% h_r = [];

lb_a = [-inf_const; -inf_const; -inf_const; -inf_const];
lb_b = [-inf_const; -inf_const; -inf_const; -inf_const];

% lb_c = x0_c;
ub_a = [inf_const; inf_const; inf_const; inf_const];
ub_b = [inf_const; inf_const; inf_const; inf_const];

vars_a = X_a{1};
vars_b = X_b{1};

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
    
    vars_a = [vars_a; U_a{t}; X_a{t+1}];
    vars_b = [vars_b; U_b{t}; X_b{t+1}];
     
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
    
    % Dynamic constraints.                
    h_a = [h_a; X_a_pred - X_a{t+1}];
    h_b = [h_b; X_b_pred - X_b{t+1}];
        
   
    % Primal variable bounds.
    lb_a = [lb_a; -inf_const; -inf_const; -inf_const; -inf_const; -inf_const; -inf_const];
    ub_a = [ub_a;  inf_const; inf_const;  inf_const; inf_const; inf_const; inf_const];
    lb_b = [lb_b; -inf_const; -inf_const; -inf_const; -inf_const; -inf_const; -inf_const];
    ub_b = [ub_b;  inf_const; inf_const; inf_const; inf_const; inf_const; inf_const];
    
    % Cost function design.
    f_a = f_a + ...
        1*(U_a{t}'*U_a{t}) + ... % Minimize acceleration
        1*((X_a{t+1}(1:2))'*(X_a{t+1}(1:2))) + 1*(X_a{t+1}(1:4)'*X_b{t+1}(1:4));
    f_b = f_b + ...
        1*(U_b{t}'*U_b{t}) + ...
         2*((X_b{t+1}(1:4) - X_a{t+1}(1:4))' * (X_b{t+1}(1:4) - X_a{t+1}(1:4)));

end

lam_a = MX.sym('lam_a', size(h_a,1));
lam_b = MX.sym('lam_b', size(h_b,1));

fudge_a = MX.sym('fudge_a', size(h_a,1));
fudge_b = MX.sym('fudge_a', size(h_b,1));

shlam_a = MX.sym('shlam_a', size(h_a,1));
shlam_b = MX.sym('shlam_b', size(h_b,1));

Lag_a = f_a; 
Lag_b = f_b;

if size(h_a,1) > 0
    Lag_a = Lag_a - lam_a'*h_a;
end

if size(h_b,1) > 0
    Lag_b = Lag_b - lam_b'*h_b;
end

gLaU = gradient(Lag_a, q_a);
gLaX = gradient(Lag_a, p_a);
gLbU = gradient(Lag_b, q_b);
gLbX = gradient(Lag_b, p_b);

% gLa = gradient(Lag_a, vars_a);
% gLb = gradient(Lag_b, vars_b);

% traj_size = size(q_a,1) + size(q_b,1) + size(p_a,1) + size(p_b,1);
% all_vars = [vars_a; vars_b; lam_a; lam_b; mu_a; mu_b];
% R = jacobian([gLaU; gLbU], [q_a; q_b]);
% Q = jacobian([gLaX; gLbX], [p_a; p_b]);
% BBAA = jacobian([h_a;h_b], [q_a;q_b;p_a;p_b]);
% B = BBAA(:,1:size(R,1));
% A = BBAA(:,size(R,1)+1:end);
% 
% fudge_constraints = [fudge_a+1; fudge_b+1];

gL = [gLaU; gLbU; gLaX; gLbX; h_a; h_b];
% mm = size(q_a,1);
% % mm = size(gLaU,1);
% nn = size(p_a,1);
% ll = size(lam_a,1);
% tt = size(vars_a,1) + size(vars_b,1);

all_vars = [q_a; q_b; p_a; p_b; lam_a; lam_b];
HH = jacobian(gL, all_vars);
% HHi = inv(HH);
% d_ua_ha = HHi(1:mm,tt+1:tt+ll);
% d_ua_hb = HHi(1:mm,tt+ll+1:tt+2*ll);
% d_ub_ha = HHi(mm+1:2*mm,tt+1:tt+ll);
% d_ub_hb = HHi(mm+1:2*mm,tt+ll+1:tt+2*ll);
% d_xa_ha = HHi(2*mm+1:2*mm+nn,tt+1:tt+ll);
% d_xa_hb = HHi(2*mm+1:2*mm+nn,tt+ll+1:tt+2*ll);
% d_xb_ha = HHi(2*mm+nn+1:2*mm+2*nn,tt+1:tt+ll);
% d_xb_hb = HHi(2*mm+nn+1:2*mm+2*nn,tt+ll+1:tt+2*ll);
% 
% d_ua_ub_xa_xb_d_ha_hb = HHi(1:tt,tt+1:tt+2*ll);
% look = HHi(1:tt,tt+2*ll+1:end);
% 
% Mat1 = (A'+Q'*inv(A)*B*inv(R)*B');
% Mat2 = [Q'*inv(A)*B*inv(R), eye(size(A,1))];
% 
% 
% Diff = d_ua_ub_xa_xb_d_ha_hb' - inv(Mat1)*Mat2;
% 
% % nn = size(p_a,1);
% % mm = size(q_a,1);
% 
% grad_fa_ub = gradient(f_a, q_b);
% grad_fa_xb = gradient(f_a, p_b);
% grad_fb_ua = gradient(f_b, q_a);
% grad_fb_xa = gradient(f_b, p_a);
% 
% grad_fb_ub = gradient(f_b, q_b);
% grad_fb_xb = gradient(f_b, p_b);
% grad_fa_ua = gradient(f_a, q_a);
% grad_fa_xa = gradient(f_a, p_a);
% 
% grads = [grad_fa_ua;
%          grad_fb_ub;
%          grad_fa_xa;
%          grad_fb_xb];
% %      
% % nn = size(A,1);
% % h_shadow = [shlam_a; shlam_b];
% 
% % h_shadow_b = shlam_b - [d_ua_hb; d_xa_hb]'*[grad_fb_ua; grad_fb_xa];
% % h_shadow_a = shlam_a - [d_ub_ha; d_xb_ha]'*[grad_fa_ub; grad_fa_xb];
% 
% h_shadow_b = shlam_b - [d_ub_hb; d_xb_hb]'*[grad_fb_ub; grad_fb_xb];
% h_shadow_a = shlam_a - [d_ua_ha; d_xa_ha]'*[grad_fa_ua; grad_fa_xa];
% % h_shadow = [shlam_a; shlam_b] - d_ua_ub_xa_xb_d_ha_hb' * grads;
% 
% % h_shadow = [shlam_a; shlam_b] - inv(A' + Q*inv(A)*B*inv(R)*B')*[Q*inv(A)*B*inv(R), eye(nn)]*grads;
% % ll = size(h_shadow,1) / 2;
% % h_shadow_a = h_shadow(1:ll,1);
% % h_shadow_b = h_shadow(ll+1:end,1);
% 
% gL = [gLaU; gLbU; gLaX; gLbX; h_a; h_b; h_shadow_a; h_shadow_b];
% all_vars = [vars_a; vars_b; lam_a; lam_b; shlam_a; shlam_b];
% HH = jacobian(gL, all_vars);




spy(HH)

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});
evalFa = Function('obj_a', {all_vars}, {f_a});
evalFb = Function('obj_b', {all_vars}, {f_b});
evalha = Function('const_a', {all_vars}, {h_a});
evalhb = Function('const_b', {all_vars}, {h_b});


% evalGrads = Function('grads', {all_vars}, {d_ua_ub_xa_xb_d_ha_hb});
% evalDiff = Function('diffs', {all_vars}, {Diff});
% evalLook = Function('look', {all_vars1}, {look});
% % evalR = Function('gamehess', {all_vars}, {R});
% % evalQ = Function('gamehess', {all_vars}, {Q});
% evalA = Function('AAA', {all_vars}, {A});
% evalB = Function('gamehess', {all_vars}, {B});
% evalhr = Function('const_r', {all_vars}, {h_r});

% evalga = Function('const_a', {all_vars}, {g_a});
% evalgb = Function('const_b', {all_vars}, {g_b});
% evalgr = Function('const_r', {all_vars}, {g_r});

extract_pa = Function('extract_pa', {all_vars}, {p_a});
extract_pb = Function('extract_pb', {all_vars}, {p_b});
% extract_pr = Function('extract_pr', {all_vars}, {vars_r});

extract_qa = Function('extract_qa', {all_vars}, {q_a});
extract_qb = Function('extract_qb', {all_vars}, {q_b});

% extract_lam_a = Function('extract_lam_a', {all_vars}, {lam_a});
% extract_lam_b = Function('extract_lam_b', {all_vars}, {lam_b});
% extract_shlam_a = Function('extract_shlam_a', {all_vars}, {shlam_a});
% extract_shlam_b = Function('extract_shlam_b', {all_vars}, {shlam_b});

% extract_mu_a = Function('extract_mu_a', {all_vars}, {mu_a});
% extract_mu_b = Function('extract_mu_b', {all_vars}, {mu_b});


lo = -inf_const*ones(size(all_vars));
hi = inf_const*ones(size(all_vars));

z0 = zeros(size(all_vars));

n = size(z0,1);
m = n;
% m = size([vars_a; vars_b;lam_a; lam_b; shlam_b; shlam_a],1);
nnz = nnz(evalH(z0));
init_loc_a = 0;
init_loc_b = size(vars_a,1);
dim_state = size(X_a{1},1);
dim_control = size(U_a{1},1);
specs = [n;m;nnz;init_loc_a;init_loc_b;dim_state;dim_control;T;num_iterations];


writez(lo, '/Users/forrest/code/packages/pathlib-master/examples/C/lo.txt');
writez(hi, '/Users/forrest/code/packages/pathlib-master/examples/C/hi.txt');
writez(z0, '/Users/forrest/code/packages/pathlib-master/examples/C/z0.txt');
write_nums(specs, '/Users/forrest/code/packages/pathlib-master/examples/C/spec.txt');

% size([vars_a; vars_b; vars_c; vars_r; lam_a; lam_b; lam_c; lam_r])
% evalH(z0)

opts = struct('with_header', true, 'casadi_int', 'int');

% C = CodeGenerator('gen.c', opts);
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

%
% system('source ~/.bash_profile');
% system('/Users/forrest/code/packages/pathlib-master/examples/C/game_shared');
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


