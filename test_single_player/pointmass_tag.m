clear; clc; close all;
import casadi.*

inf_const = 999999999999.9;

dt = 0.25;
T = 30;
num_iterations = 80;

constraint_radius = 2;

X_a{1} = MX.sym(['X_a_' num2str(0)], 4);
X_b{1} = MX.sym(['X_b_' num2str(0)], 4);

a_bound_a = 1; % m/s^2
a_bound_b = 1;
v_bound = 3;

x0_a = [-2.5;1.1;-.1;-.1];
x0_b = [8.1;5;0;0];

f_a = 0;
f_b = 0;

h_a = [];
h_b = [];


lb_a = x0_a;
lb_b = x0_b;

ub_a = x0_a;
ub_b = x0_b;

g_a = [];
g_b = [];

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
    lb_a = [lb_a; -a_bound_a; -a_bound_a; -inf_const; -inf_const; -v_bound; -v_bound];
    ub_a = [ub_a;  a_bound_a; a_bound_a;  inf_const; inf_const; v_bound; v_bound];
    lb_b = [lb_b; -a_bound_b; -a_bound_b; -inf_const; -inf_const; -v_bound; -v_bound];
    ub_b = [ub_b;  a_bound_b; a_bound_b; inf_const; inf_const; v_bound; v_bound];
    
    g_a = [g_a; 
        (X_a{t+1}(1:2))'*(X_a{t+1}(1:2)) - constraint_radius^2;
        (X_a{t+1}(1:2)-X_b{t+1}(1:2))'*(X_a{t+1}(1:2)-X_b{t+1}(1:2)) - 1];

    g_b = [g_b; 
        (X_b{t+1}(1:2))'*(X_b{t+1}(1:2)) - constraint_radius^2];
    % Cost function design.
    f_a = f_a + ...
        1*(U_a{t}(1:2)'*U_a{t}(1:2)) + ... % Minimize acceleration
        -.001*(X_a{t+1}(1:2) - X_b{t+1}(1:2))' * (X_a{t+1}(1:2) - X_b{t+1}(1:2)) + ...
        .1*(X_a{t+1}(1:2)'*X_a{t+1}(1:2));
    f_b = f_b + ...
        .1*(U_b{t}(1:2)'*U_b{t}(1:2)) + ...
         1*(X_b{t+1}(1:2) - X_a{t+1}(1:2))' * (X_b{t+1}(1:2) - X_a{t+1}(1:2));

end

lam_a = MX.sym('lam_a', size(h_a,1));
lam_b = MX.sym('lam_b', size(h_b,1));
% lam_r = MX.sym('lam_r', size(h_r,1));

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

% if size(h_r,1) > 0
%     Lag_r = Lag_r - lam_r'*h_r;
% end
% if size(g_r,1) > 0
%     Lag_r = Lag_r - mu_r'*g_r;
% %     Lag_a = Lag_a - mu_ra'*g_r;
% %     Lag_b = Lag_b - mu_rb'*g_r;
% %     Lag_c = Lag_b - mu_rc'*g_r;
% %     mu_a = [mu_a;mu_ra];
% %     mu_b = [mu_b;mu_rb];
% %     mu_c = [mu_c;mu_rc]; 
% end



gLa = gradient(Lag_a, vars_a);
gLb = gradient(Lag_b, vars_b);
% gLr = gradient(Lag_r, vars_r);

gL = [gLa; gLb; h_a; h_b; g_a; g_b];

all_vars = [vars_a; vars_b; lam_a; lam_b; mu_a; mu_b];
HH = jacobian(gL, all_vars);
% spy(HH)

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});
evalFa = Function('obj_a', {all_vars}, {f_a});
evalFb = Function('obj_b', {all_vars}, {f_b});
evalha = Function('const_a', {all_vars}, {h_a});
evalhb = Function('const_b', {all_vars}, {h_b});
% evalhr = Function('const_r', {all_vars}, {h_r});

evalga = Function('const_a', {all_vars}, {g_a});
evalgb = Function('const_b', {all_vars}, {g_b});
% evalgr = Function('const_r', {all_vars}, {g_r});

extract_pa = Function('extract_pa', {all_vars}, {p_a});
extract_pb = Function('extract_pb', {all_vars}, {p_b});
% extract_pr = Function('extract_pr', {all_vars}, {vars_r});

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

% ones_r = ones(size(vars_r,1),1);
% ones_lam_r = ones(size(lam_r,1),1);
% ones_mu_r = ones(size(mu_r,1),1);

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
    % Stay in lane params.
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

% figure; hold on;
% plot([-2.5,-2.5],[-5,500],'k');
% plot([2.5,2.5],[-5,500],'k');
% plot([7.5,7.5],[-5,500],'k');
% 
% for t = 1:T+1
%     hold on;
%     plot(YA(t,1),XA(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
%     [Xsrc_a, Ysrc_a] = drawShearedRecircle(XA(t,1),YA(t,1),shear,half_w,half_h);
%     ca = plot(Xsrc_a,Ysrc_a,'b');
% %     ca1 = viscircles([YA(t,1),XA(t,1)]+off_for,bounding_radius,'Color', 'b');
% %     ca2 = viscircles([YA(t,1),XA(t,1)]+off_bac,bounding_radius,'Color', 'b');
%     plot(YB(t,1),XB(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%     [Xsrc_b, Ysrc_b] = drawShearedRecircle(XB(t,1),YB(t,1),shear,half_w,half_h);
%     cb = plot(Xsrc_b,Ysrc_b,'r');
% %     cb1 = viscircles([YB(t,1),XB(t,1)]+off_for,bounding_radius,'Color', 'r');
% %     cb2 = viscircles([YB(t,1),XB(t,1)]+off_bac,bounding_radius,'Color', 'r');
%     plot(YC(t,1),XC(t,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
%     [Xsrc_c, Ysrc_c] = drawShearedRecircle(XC(t,1),YC(t,1),shear,half_w,half_h);
%     cc = plot(Xsrc_c,Ysrc_c,'g');
% %     cc1 = viscircles([YC(t,1),XC(t,1)]+off_for,bounding_radius,'Color', 'g');
% %     cc2 = viscircles([YC(t,1),XC(t,1)]+off_bac,bounding_radius,'Color', 'g');
%     pos = XB(t,1);
%     axis([-30 30 pos-30 pos+30]);
%     pause(0.2);
%     delete(ca);
%     delete(cb);
%     delete(cc);
% %     delete(ca2);
% %     delete(cb2);
% %     delete(cc2);
% end


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
[status,cmdout] = system(cmd);

%
% system('source ~/.bash_profile');
% system('/Users/forrest/code/packages/pathlib-master/examples/C/game_shared');
%%
ZZ = readz('/Users/forrest/code/packages/pathlib-master/examples/C/z.csv');

writerObj = VideoWriter('pointmass_tag','MPEG-4');
writerObj.FrameRate = 8;
open(writerObj);
% Z = peaks; surf(Z); 
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
hold on;
% for k = 1:20 
% surf(sin(2*pi*k/20)*Z,Z)
% frame = getframe;
% writeVideo(writerObj,frame);
% end
% close(writerObj);

% figure; hold on;


viscircles([0,0],constraint_radius,'color','k');
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

    % P_AB_X = PR(1:12:end); P_AB_Y = PR(2:12:end);
    % P_BA_X = PR(3:12:end); P_BA_Y = PR(4:12:end);
    % P_AC_X = PR(5:12:end); P_AC_Y = PR(6:12:end);
    % P_CA_X = PR(7:12:end); P_CA_Y = PR(8:12:end);
    % P_BC_X = PR(9:12:end); P_BC_Y = PR(10:12:end);
    % P_CB_X = PR(11:12:end); P_CB_Y = PR(12:12:end);

    % plot([-2.5,-2.5],[-5,500],'k');
    % plot([2.5,2.5],[-5,500],'k');
    % plot([7.5,7.5],[-5,500],'k');
   
    path_a = plot(XA(:,1), YA(:,1),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
    path_b = plot(XB(:,1),YB(:,1),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'r');

    spot_a = plot(XA(1,1), YA(1,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    spot_b = plot(XB(1,1),YB(1,1),'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     cur_a = viscircles([XA(1,1),YA(1,1)],1,'Color', 'b');
%     cur_b = viscircles([XB(1,1),YB(1,1)],1,'Color', 'r');
%     time = text(XA(1,1)+15, YA(1,1)+15,string(iter*dt));
    axis([-15 15 -15 15]);
    frame = getframe;
    writeVideo(writerObj,frame);
    if num_iterations > 1
%         pause(dt);
        delete(path_a);
        delete(path_b);
        delete(spot_a);
        delete(spot_b);
%         delete(time);
    end
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

close(writerObj);
