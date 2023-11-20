clear; clc; close all;
import casadi.*

inf_const = 999999999999.9;

dt = 0.025;
T = 200;
num_iterations = 1;

X{1} = MX.sym(['X_' num2str(1)], 6);


mu = 0.1;
u_bound = 10;
l = 4;
K = 100;
C = 5;

g = 7;
h = 10;
x0 = [0;h;h-l;0;-15;-15]; % x, yc, yp, vx, vyc, vyp


f_x = 0;

h_x = [X{1}-x0];
g_x = [];
vars_x = X{1};

states = X{1};
controls = [];
lams = [];
gams = [];

for t = 1:T
    % Introduce symbolic vars.
    U{t} = MX.sym(['U_' num2str(t)], 1);
    Y{t} = MX.sym(['Y_' num2str(t)], 1);
    Z{t} = MX.sym(['Z_' num2str(t)], 1);
    X{t+1} = MX.sym(['X_' num2str(t+1)], 6);

    vars_x = [vars_x; U{t}; X{t+1}];
    vars_y{t} = [Y{t}; X{t+1}];
    vars_z{t} = [Z{t}; X{t+1}];
   
    % Place vars in lists for easy access later.
    states = [states; X{t+1}];
    controls = [controls; U{t}];
    lams = [lams; Y{t}];
    gams = [gams; Z{t}];
    
    % Dynamic equations
    X_pred = X{t} + dt * [X{t}(4); 
                          X{t}(5); 
                          X{t}(6); 
                          U{t}(1)+Z{t};
                          -K*(X{t}(2)-X{t}(3)-l)-C*(X{t}(5)-X{t}(6))-g;
                           K*(X{t}(2)-X{t}(3)-l)+C*(X{t}(5)-X{t}(6))+Y{t}-g];
    
    % Dynamic constraints.                
    h_x = [h_x; X_pred - X{t+1}];
    h_y{t} = [X_pred - X{t+1}];
    h_z{t} = [X_pred - X{t+1}];
    if (t>1)
        h_y{t-1} = [h_y{t-1}; X_pred-X{t+1}];
%         h_y{t-1} = [h_y{t-1}; X_pred-X{t+1}];
    end


%     g_x = g_x;
    g_y{t} = X{t+1}(3);
    g_z{t} = 0.00000001+Y{t}*Y{t}*mu^2 - Z{t}*Z{t};
%     g_y{t} = [];
%     g_z{t} = [];

    f_x = f_x + 1.0*(U{t}'*U{t}) + (X{t}(1)-5)^2;
    
    f_y{t} = Y{t}*Y{t};
    
    f_z{t} = X{t+1}(4)^2;

end

for t = 1:T

    if (t < T)
        g_y{t} = X{t+2}(3);
        vars_y{t} = [vars_y{t}; X{t+2}];
%         vars_z{t} = [vars_z{t}; X{t+2}];
    else
        g_y{t} = [];
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
%     gLy{t} = gradient(Lag_y{t}, vars_y{t}); % n + l
%     gLz{t} = gradient(Lag_z{t}, vars_z{t}); % n + g
end


gL = [gLx; all_gLy; all_gLz];
% for t = 1:T
%     gL = [gL; gLy{t}; gLz{t}];
% end
gL = [gL; h_x; g_x];
all_g_y = [];
all_g_z = [];
for t = 1:T
    all_g_y = [all_g_y; g_y{t}];
    all_g_z = [all_g_z; g_z{t}];
%     gL = [gL; g_y{t}; g_z{t}];
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
%     all_vars = [all_vars; lam_y{t}; lam_z{t}];
end
all_vars = [all_vars; all_lam_y; all_lam_z];
all_vars = [all_vars; mu_x];
for t = 1:T
    all_mu_y = [all_mu_y; mu_y{t}];
    all_mu_z = [all_mu_z; mu_z{t}];
%     all_vars = [all_vars; mu_y{t}; mu_z{t}];
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


lo = [-inf_const*ones(size(states));
    -u_bound*ones(size(controls));
    -inf_const*ones(size(lams));
    -inf_const*ones(size(gams));
    -inf_const*ones(size(lam_x));
    -inf_const*ones(size(all_lam_y));
    -inf_const*ones(size(all_lam_z));
    zeros(size(mu_x));
    zeros(size(all_mu_y));
    zeros(size(all_mu_z))];
hi = [inf_const*ones(size(states));
    u_bound*ones(size(controls));
    inf_const*ones(size(lams));
    inf_const*ones(size(gams));
    inf_const*ones(size(lam_x));
    inf_const*ones(size(all_lam_y));
    inf_const*ones(size(all_lam_z));
    inf_const*ones(size(mu_x));
    inf_const*ones(size(all_mu_y));
    inf_const*ones(size(all_mu_z))];

z0 = zeros(size(all_vars,1),1);

n = size(z0,1);
m = size([states; controls; lams; gams; lam_x; all_lam_y; all_lam_z], 1);
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
lams = full(extract_lams(z));
gams = full(extract_gams(z));

x = states(1:6:end);
yc = states(2:6:end);
yp = states(3:6:end);
vx = states(4:6:end);
vc = states(5:6:end);
vp = states(6:6:end);
ux = controls;

% writerObj = VideoWriter('bouncing_ball_game','MPEG-4');
% writerObj.FrameRate = 20;
% open(writerObj);
% axis tight
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
hold on;


plot([-10,10],[0,0],'-k')
path_a = plot([5], [0],'x', 'MarkerSize', 20, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
axis('equal')
for t=1:T+1
    pp = fimplicit(@(zx,zy) (zx-x(t))^2/(l/2)^2+(zy-0.5*(yc(t)+yp(t)))^2/((yc(t)-yp(t))/2)^2-1, '-b');
    axis([-10 10 -5 +15]);
%     pp = ezplot('(x-x(t)+(y-0.5*(yc(t)+yp(t)))^2/((yc(t)-yp(t))/2)^2=1');
%     pc = plot(x(t),yc(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
%     pp = plot(x(t),yp(t),'or', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'g');
%     frame = getframe;
%     writeVideo(writerObj,frame);
    pause(dt);
%     delete(pc);
    delete(pp);
end

% close(writerObj);
hold on;




