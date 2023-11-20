clear; clc; close all;
import casadi.*

inf_const = 999999999999999;

dt = 0.04;
T = 30;
num_iterations = 150;

X{1} = MX.sym(['X_' num2str(1)], 8);

wall_angle = pi/4;

Mc = 1;
Me = 0.1;
bar = 6;

mu = 0.1;
u_bound = inf_const;
l = 4;
K = 100;
C = 1;

g = 9.81;
h0 = 25;

x0 = [-10; %cx
    h0;   %cy
    -1;   %vx
    -3;   %vy
    l;   %dl
    0;   %vl
    l;   %dr
    0];  %vr

f_x = 0;

h_x = [];
g_x = [];
vars_x = X{1};

states = X{1};
controls = [];
lams = [];
gams = [];

for t = 1:T
    % Introduce symbolic vars.
    U{t} = MX.sym(['U_' num2str(t)], 2);
    Y{t} = MX.sym(['normal_force_' num2str(t)], 2);
    Z{t} = MX.sym(['tangential_force_' num2str(t)], 2);
    
    X{t+1} = MX.sym(['X_' num2str(t+1)], 8);

    vars_x = [vars_x; U{t}; X{t+1}];
    vars_y{t} = [Y{t}; X{t+1}];
    vars_z{t} = [Z{t}; X{t+1}];
   
    % Place vars in lists for easy access later.
    states = [states; X{t+1}];
    controls = [controls; U{t}];
    lams = [lams; Y{t}];
    gams = [gams; Z{t}];
                            
    left_spring_force = K*(X{t}(5)-l) + C*X{t}(6);
    right_spring_force = K*(X{t}(7)-l) + C*X{t}(8);
    
    left_com_accel = -[sin(wall_angle); cos(wall_angle)] * left_spring_force/(Mc+Me) + [cos(wall_angle); -sin(wall_angle)] * Z{t}(1) / (Mc+2*Me);
    right_com_accel = -[-sin(wall_angle); cos(wall_angle)] * right_spring_force/(Mc+Me) + [cos(wall_angle); sin(wall_angle)] * Z{t}(2) / (Mc+2*Me);
    
    left_length_accel = (-Y{t}(1) - left_spring_force) / Me;
    right_length_accel = (-Y{t}(2) - right_spring_force) / Me;
    
    X_pred = X{t} + dt * [X{t}(3)+dt*(left_com_accel(1)+right_com_accel(1)+U{t}(1));
                          X{t}(4)+dt*(left_com_accel(2)+right_com_accel(2)-g+U{t}(2));
                          left_com_accel+right_com_accel+[U{t}(1);-g+U{t}(2)];
                          X{t}(6)+dt*left_length_accel;
                          left_length_accel;
                          X{t}(8)+dt*right_length_accel;
                          right_length_accel];

    % Dynamic constraints.                
    h_x = [h_x; X_pred - X{t+1}];
    h_y{t} = [X_pred - X{t+1}];
    h_z{t} = [X_pred - X{t+1}];
    
    point_left = [X{t+1}(1); X{t+1}(2)] + X{t+1}(5) * [-sin(wall_angle); -cos(wall_angle)];
    point_right = [X{t+1}(1); X{t+1}(2)] + X{t+1}(7) * [sin(wall_angle); -cos(wall_angle)];
    
    if (t>1)
        g_x = [g_x; X{t+1}(2) - l - bar];
    end
    
    g_y{t} = [point_left(2) - tan(-wall_angle)*point_left(1);
              point_right(2) - tan(wall_angle)*point_right(1)];
    
%     g_z{t} = [(0.01)^2+Y{t}(1)^2*mu^2 - Z{t}(1)^2;
%               (0.01)^2+Y{t}(2)^2*mu^2 - Z{t}(2)^2];
          
    g_z{t} = [0.01+mu*Y{t} - Z{t};
              0.01+mu*Y{t} + Z{t}];

          
    if (t+1>=T)
        const = 0;
    else
        const = 0;
    end
    f_x = f_x + 10.0*(U{t}'*U{t}) - const*(X{t+1}(2)) - 50*abs(X{t+1}(3));

           
    f_y{t} = Y{t}'*Y{t};
    
    f_z{t} = ([X{t+1}(3);X{t+1}(4)]'*[cos(wall_angle); sin(wall_angle)])^2 + ...
             ([X{t+1}(3);X{t+1}(4)]'*[cos(wall_angle); -sin(wall_angle)])^2;

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
all_gLy = [];
all_gLz = [];
for t = 1:T
    lam_y{t} = MX.sym('lam_y', size(h_y{t},1));
    lam_z{t} = MX.sym('lam_z', size(h_z{t},1));
    mu_y{t} = MX.sym('mu_y', size(g_y{t},1));
    mu_z{t} = MX.sym('mu_z', size(g_z{t},1));

    Lag_y{t} = f_y{t};
    Lag_z{t} = f_z{t};


    if size(h_y{t},1) > 0
        Lag_y{t} = Lag_y{t} - lam_y{t}'*h_y{t};
    end
    if size(g_y{t},1) > 0
        Lag_y{t} = Lag_y{t} - mu_y{t}'*g_y{t};
    end
    if size(h_z{t},1) > 0
        Lag_z{t} = Lag_z{t} - lam_z{t}'*h_z{t};
    end
    if size(g_z{t},1) > 0
        Lag_z{t} = Lag_z{t} - mu_z{t}'*g_z{t};
    end
    all_gLy = [all_gLy; gradient(Lag_y{t}, vars_y{t})];
    all_gLz = [all_gLz; gradient(Lag_z{t}, vars_z{t})];
end


gL = [gLx; all_gLy; all_gLz];
gL = [gL; h_x; g_x];
all_g_y = [];
all_g_z = [];
for t = 1:T
    all_g_y = [all_g_y; g_y{t}];
    all_g_z = [all_g_z; g_z{t}];
end
gL = [gL; all_g_y; all_g_z];

all_vars = [states; controls; lams; gams; lam_x];
all_lam_y = [];
all_lam_z = [];
all_mu_y = [];
all_mu_z = [];
for t = 1:T
    all_lam_y = [all_lam_y; lam_y{t}];
    all_lam_z = [all_lam_z; lam_z{t}];
end
all_vars = [all_vars; all_lam_y; all_lam_z];
all_vars = [all_vars; mu_x];
for t = 1:T
    all_mu_y = [all_mu_y; mu_y{t}];
    all_mu_z = [all_mu_z; mu_z{t}];
end
all_vars = [all_vars; all_mu_y; all_mu_z];

HH = jacobian(gL, all_vars);
% spy(HH)

evalH = Function('gamehess', {all_vars}, {HH});
evalG = Function('gamegrad', {all_vars}, {gL});

extract_states = Function('extract_states', {all_vars}, {states});
extract_controls = Function('extract_controls', {all_vars}, {controls});
extract_lams = Function('extract_lams', {all_vars}, {lams});
extract_gams = Function('extract_gams', {all_vars}, {gams});


lo = [x0;
    -inf_const*ones(size(states,1)-8,1);
    -u_bound*ones(size(controls));
    -inf_const*ones(size(lams));
    -inf_const*ones(size(gams));
    -inf_const*ones(size(lam_x));
    -inf_const*ones(size(all_lam_y));
    -inf_const*ones(size(all_lam_z));
    zeros(size(mu_x));
    zeros(size(all_mu_y));
    zeros(size(all_mu_z))];
hi = [x0;
    inf_const*ones(size(states,1)-8,1);
    u_bound*ones(size(controls));
    inf_const*ones(size(lams));
    inf_const*ones(size(gams));
    inf_const*ones(size(lam_x));
    inf_const*ones(size(all_lam_y));
    inf_const*ones(size(all_lam_z));
    inf_const*ones(size(mu_x));
    inf_const*ones(size(all_mu_y));
    inf_const*ones(size(all_mu_z))];

z0 = zeros(size(all_vars,1), 1);
z0(1:8) = x0;

n = size(z0,1);
m = size([states; controls; lams; gams; lam_x; all_lam_y; all_lam_z], 1);
nnz = nnz(evalH(z0));
init_loc_a = size(states,1);
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
disp 'Generated game binaries.';
str = ['Simulating ' num2str(num_iterations*dt) ' seconds of game play.'];
disp(str);

cmd = "cd /Users/forrest/code/packages/pathlib-master/examples/C && " + ...
      "export DYLD_LIBRARY_PATH=/Users/forrest/code/packages/pathlib-master/lib/osx && " + ...
      "export PATH_LICENSE_STRING='2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0' && " + ...
      "cp ~/Documents/MATLAB/fb_games/test_single_player/gen.c . && make && ./game_shared";
[status,cmdout] = system(cmd)

%%
ZZ = readz('/Users/forrest/code/packages/pathlib-master/examples/C/z.csv');

writerObj = VideoWriter('wallball','MPEG-4');
writerObj.FrameRate = 20;
open(writerObj);
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
hold on;

% figure;
% hold on;

bound = 15;
bound = bound + l + 3;



plot([-bound,bound],[bar, bar], '--r');

h=fill([-bound -bound bound bound],[-bound bar bar -bound],'red');
h.FaceAlpha=0.5;
plot([-bound,0],tan(-wall_angle)*[-bound,0],'-k','linewidth', 4)
plot([0,bound],tan(wall_angle)*[0,bound],'-k','linewidth', 4)
% path_a = plot([5], [0],'x', 'MarkerSize', 20, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
axis('equal')

mpc_on = (num_iterations > 1);

for iter = 1:num_iterations
    z = ZZ((iter-1)*n+1:iter*n);
    states = full(extract_states(z));
    controls = full(extract_controls(z));
    lams = full(extract_lams(z));
    gams = full(extract_gams(z));

    x = states(1:8:end);
    y = states(2:8:end);
    vx = states(3:8:end);
    vy = states(4:8:end);
    dl = states(5:8:end);
    vl = states(6:8:end);
    dr = states(7:8:end);
    vr = states(8:8:end);

    normal_left = lams(1:2:end);
    normal_right = lams(2:2:end);
    friction_left = gams(1:2:end);
    friction_right = gams(2:2:end);

    ux = controls;
    control_right = ux(1:2:end);
    control_up = ux(2:2:end);

    if (mpc_on)
        t = 1;
        pp = fimplicit(@(zx,zy) ((zx-x(t))*cos(-wall_angle) + (zy-y(t))*sin(-wall_angle))^2 / (dr(t)^2) + ((zx-x(t))*sin(-wall_angle) - (zy-y(t))*cos(-wall_angle))^2 / (dl(t)^2) - 1, '-b');
        point_center = [x(t); y(t)];
        v_point_left(t,:) = point_center + dl(t) * [-sin(wall_angle); -cos(wall_angle)];
        v_point_right(t,:) = point_center + dr(t) * [sin(wall_angle); -cos(wall_angle)];

        arrow_up = quiver(x(t), y(t)-control_up(t),0, control_up(t),0,'-g','linewidth',2);
        arrow_right = quiver(x(t)-control_right(t), y(t),control_right(t), 0,0,'-g','linewidth',2);

       
        traj = plot(x, y,'or', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
        axis([-bound bound -0.5*bound 1.5*bound]);
        frame = getframe;
        writeVideo(writerObj,frame);
%         pause(dt);
        delete(pp);
        delete(arrow_up);
        delete(arrow_right);
        delete(traj);
    else
        for t=1:T+1
            pp = fimplicit(@(zx,zy) ((zx-x(t))*cos(-wall_angle) + (zy-y(t))*sin(-wall_angle))^2 / (dr(t)^2) + ((zx-x(t))*sin(-wall_angle) - (zy-y(t))*cos(-wall_angle))^2 / (dl(t)^2) - 1, '-b');
            point_center = [x(t); y(t)];
            v_point_left(t,:) = point_center + dl(t) * [-sin(wall_angle); -cos(wall_angle)];
            v_point_right(t,:) = point_center + dr(t) * [sin(wall_angle); -cos(wall_angle)];
            if (t<T+1)
                arrow_up = quiver(x(t), y(t)-control_up(t),0, control_up(t),0,'-b','linewidth',2);
                arrow_right = quiver(x(t)-control_right(t), y(t),control_right(t), 0,0,'-b','linewidth',2);
            end

            px = [v_point_left(t,1) point_center(1) v_point_right(t,1)];
            py = [v_point_left(t,2) point_center(2) v_point_right(t,2)];
            axis([-bound bound -0.5*bound 1.5*bound]);
        %     frame = getframe;
        %     writeVideo(writerObj,frame);
            pause(dt);
            delete(pp);
            if (t < T+1)
                delete(arrow_up);
                delete(arrow_right);
            end
        end
    end
end

close(writerObj);
% hold on;




