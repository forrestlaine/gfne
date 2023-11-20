%% Is there a car?
clear; clc; close all;

global dt;
dt = 0.2;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;
global lane_width;
global draw_borders;
draw_borders = true;


lane_width = 4;
border_left = -lane_width;
border_right = lane_width;

global speed_pen lat_speed_pen;

global R;


global half_l half_w;
global buffer_front;
global buffer_front_tip;
global buffer_behind_start;
global buffer_behind_length_a;
global buffer_behind_length_b;
global buffer_behind_length_c;
global buffer_width_wide;
global buffer_width_narrow;
half_l = 1.5;
half_w = 0.75;



buffer_front = 2;

buffer_behind_start = 18;
buffer_behind_length_a = 1.3;
buffer_behind_length_b = 1.5;
buffer_behind_length_c = 1;
buffer_front_tip = 3;
buffer_width_wide = 1.8;
buffer_width_narrow = 1.7;


R = diag([.001,.001]);
speed_pen = 2;
lat_speed_pen = .75;

global gv_a1 gv_b1 gv_c1 gv_d1;



global T
T = 50;
mpc_iters = 200;

gv_a1 = 25;
gv_b1 = 27;
gv_c1 = 23;
% gv_d1 = 23;

x0 = [border_right/2; 0;0;  gv_a1;
      border_left/2; -18;0; gv_b1;%1
%       border_left/2; 18;0;  gv_b1;%2
      border_right/2; 25;0; gv_c1];%3
  
% x0 = [1.7479;
%   130.7328;
%    -0.0381;
%    23.3081;
%    -2.0000;
%   130.9951;
%          0;
%    27.0000;
%    -2.0000;
%    90.9951;
%          0;
%    27.0000;
%     2.0000;
%   154.0049;
%          0;
%    23.0000];


rollout_controls{1} = [0.01; -.1];
rollout_controls{2} = [0; 0];
rollout_controls{3} = [0; 0];
% rollout_controls{4} = [0; 0];
  
m{1} = 2;
m{2} = 2;
m{3} = 2;
% m{4} = 2;

full_scope{1} = true;
full_scope{2} = true;
full_scope{3} = true;
% full_scope{4} = true;

N = 3;
n = N*4;


for t = 1:T
    
    lf{t,1} = @running_cost_a1;
    lf{t,2} = @running_cost_b1;
    lf{t,3} = @running_cost_c1;
%     lf{t,4} = @running_cost_d1;
  
    gf{t} = assign_running_ineq_constraints();
    gi{t} = assign_running_ineq_ownerships();
    
    h{t} = @empty_constraint;
    hi{t} = {};

end
lf{T+1,1} = @empty_cost;
lf{T+1,2} = @empty_cost;
lf{T+1,3} = @empty_cost;
% lf{T+1,4} = @empty_cost;

gf{T+1} = assign_terminal_ineq_constraints();
gi{T+1} = assign_terminal_ineq_ownerships();

h{T+1} = @empty_constraint;
hi{T+1} = {};
% hi{T+1}{1} = 2;
% hi{T+1}{2} = 3;
% hi{T+1}{3} = 3;


%% Run 

params.use_ws = false;
params.open_loop = true;
params.use_initialization = false;
params.debug_plot = true;
params.linear = false;
params.full_scope = full_scope;
params.print_frequency = 10;
params.rollout_controls = rollout_controls;
params.wrap_val = 36;

[xvals,durs,iters] = mpgp_fb_game_solver_active_set_shared(@f, h, hi, gf, gi, lf, n, m, N, T, mpc_iters, x0, params, @plot_solution);
solution_set{1} = xvals;
durations{1} = durs;

%%
draw_borders = false;
S = numel(xvals);
for k = 1:S
    xval{k} = xvals{k}{1};
end

for k = S+2:S+T+1
    xval{k-1} = xvals{S}{k-S};
end

plot_solution(xval,0)
% draw_borders = true;
% plot_solution(xval,0);
%%

load('can_merge');
plot_solution(xval,0);


%% Plotting

function plot_solution(xval,x0)

    make_vid = true;
    make_frames = false;
    close all;
    fh = figure('Position', [300 300 150 750]); 
    if make_vid
        writerObj = VideoWriter('two_lane_can_merge_no_borders','MPEG-4');
        writerObj.FrameRate = 20;
        open(writerObj);
%         axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        set(gcf,'Color', 'w');%[.81,.81,.81]);
    end
    hold on;
    global T dt draw_borders buffer_behind_length_a buffer_behind_length_b buffer_behind_length_c;
    draw_borders = false;
    draw_lanes();

    for t = 1:numel(xval)
        xx(t,:) = xval{t};
        
%         spot_b1 = viscircles(xx(t,5:6),car_r,'Color', [1,.7,0],'LineStyle',':');
        if draw_borders
            border_a1 = draw_boundary(xx(t,1:2), buffer_behind_length_a);
            border_b1 = draw_boundary(xx(t,5:6), buffer_behind_length_b);
            border_c1 = draw_boundary(xx(t,9:10), buffer_behind_length_c);
%             border_d1 = draw_boundary(xx(t,13:14), buffer_behind_length_c);
        end
%         
%         color_a = [1,.7,0];
%         color_b = [1,.7,0];
%         color_c = [1,.7,0];
        color_a = 'b';
        color_b = 'g';
        color_c = [1,.7,0];
        car_a1 = draw_car(xx(t,1:2),color_a,1,pi/2);
        car_b1 = draw_car(xx(t,5:6),color_b,1,pi/2);
        car_c1 = draw_car(xx(t,9:10),color_c,1,pi/2);
%         car_d1 = draw_car(xx(t,13:14),color_c,1,pi/2);
        

                
%         spot_a1 = viscircles(xx(t,1:2),car_r,'Color', 'b','LineStyle',':');
        
        
        
%         spot_c1 = viscircles(xx(t,9:10),car_r,'Color', 'g','LineStyle',':');
%         border_c1 = draw_boundary(xx(t,9:10));
        
        
        
        
        axis([-4,+4,xx(t,2)-35,xx(t,2)+35]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        elseif make_frames
            axis([-4,+4,xx(t,6)-65,xx(t,6)+20]);
            axis('equal');
            if mod(t,floor(T/5)) == 1
                name = ['gfne_polite_frame_' num2str(t) '.png'];
                saveas(gcf,name)
            end
        else
            pause(0.1);
        end
        delete(car_a1);
        delete(car_b1);
        delete(car_c1);
%         delete(car_d1);
        if draw_borders
            delete(border_a1);
            delete(border_b1);
            delete(border_c1);
%             delete(border_d1);
        end
        
    end
    if make_vid
        close(writerObj);
    end
end

%%

function const = empty_constraint(~)
    const = zeros(0,1);
end

function val = empty_cost(~)
    val = 0;
end

%% Running Constraints

function const = no_lat(x)
    const = x(14);
end

function const = running_ineq_constraint_a1(x)
    % collision zone of b
    global buffer_behind_length_a buffer_behind_length_c;
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    const = avoid_polygon(xa1,xb1,buffer_behind_length_a,buffer_behind_length_c);
end

function const = running_ineq_constraint_a2(x)
    % collision zone of xc
    global buffer_behind_length_a buffer_behind_length_b;
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    const = avoid_polygon(xa1,xc1,buffer_behind_length_a,buffer_behind_length_b);
end

function const = running_ineq_constraint_a3(x)
    % collision zone of xd
    global buffer_behind_length_b buffer_behind_length_c;
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    xd1 = x(13:14);
    const = avoid_polygon(xa1,xd1,buffer_behind_length_c, buffer_behind_length_b);
end

function const = running_ineq_constraint_a4(x)
    global lane_width half_w;
    xa1 = x(1:2);
    const = left_lane_ineq(xa1);
end

function const = running_ineq_constraint_a5(x)
    global lane_width half_w;
    xa1 = x(1:2);
    const = right_lane_ineq(xa1);
end

function const = running_ineq_constraint_a6(x)
    global lane_width half_w;
    xb1 = x(5:6);
    const = left_lane_ineq(xb1);
end

function const = running_ineq_constraint_a7(x)
    global lane_width half_w;
    xb1 = x(5:6);
    const = right_lane_ineq(xb1);
end

function const = running_ineq_constraint_a8(x)
    global lane_width half_w;
    xc1 = x(9:10);
    const = left_lane_ineq(xc1);
end

function const = running_ineq_constraint_a9(x)
    global lane_width half_w;
    xc1 = x(9:10);
    const = right_lane_ineq(xc1);
end

function const = running_ineq_constraint_a10(x)
    global lane_width half_w;
    xd1 = x(13:14);
    const = left_lane_ineq(xd1);
end

function const = running_ineq_constraint_a11(x)
    global lane_width half_w;
    xd1 = x(13:14);
    const = right_lane_ineq(xd1);
end

function const = left_lane_ineq(x)
    global lane_width half_w;
    const = [x(1)-1.9*half_w+(lane_width)];
end

function const = right_lane_ineq(x)
    global lane_width half_w;
    const = [-x(1)-1.9*half_w+(lane_width)];
end

function const = assign_running_ineq_constraints()
    const{1} = @running_ineq_constraint_a1;
    const{2} = @running_ineq_constraint_a2; 
%     const{3} = @running_ineq_constraint_a3;
    const{3} = @running_ineq_constraint_a4; 
    const{4} = @running_ineq_constraint_a5;
%     const{6} = @running_ineq_constraint_a6; 
%     const{7} = @running_ineq_constraint_a7;
%     const{8} = @running_ineq_constraint_a8; 
%     const{9} = @running_ineq_constraint_a9;
%     const{10} = @running_ineq_constraint_a10; 
%     const{11} = @running_ineq_constraint_a11;
end

function ownerships = assign_running_ineq_ownerships()
    ownerships{1} = [1];  % front and back
    ownerships{2} = [1];  % front and middle
%     ownerships{3} = [1];  % middle and back
    ownerships{3} = 1;
    ownerships{4} = 1;
%     ownerships{6} = 2;
%     ownerships{7} = 2;
%     ownerships{8} = 3;
%     ownerships{9} = 3;
%     ownerships{10} = 4;
%     ownerships{11} = 4;
end


function const = avoid_polygon(va, vb,length_a,length_b)
    global buffer_width_wide buffer_width_narrow buffer_front buffer_behind_start buffer_front_tip;
    % va avoids vb
    
    % helper vertices
    va_tr = [buffer_width_wide; buffer_front];
    va_tip = [0; buffer_front_tip];
    va_tl = [-buffer_width_narrow; buffer_front];
    va_bl = [-buffer_width_wide; -buffer_behind_start];
    va_br = [buffer_width_narrow; -buffer_behind_start-length_a];
    va = va(1:2);
    
    vb_tr = vb(1:2)+[buffer_width_wide; buffer_front];
    vb_tip = vb(1:2)+[0; buffer_front_tip];
    vb_tl = vb(1:2)+[-buffer_width_narrow; buffer_front];
    vb_bl = vb(1:2)+[-buffer_width_wide; -buffer_behind_start];
    vb_br = vb(1:2)+[buffer_width_narrow; -buffer_behind_start-length_b];
    
    % collision region vertices
    v1 = vb_tip - va_br;
    v2 = vb_tip - va_bl;
    v3 = vb_tr - va_bl;
    v4 = vb_tr - va_tl;
    v5 = vb_br - va_tl;
    v6 = vb_br - va_tip;
    v7 = vb_bl - va_tip;
    v8 = vb_bl - va_tr;
    v9 = vb_tl - va_tr;
    v10 = vb_tl - va_br;
    
    % collision region faces
    v12 = v1-v2;
    v23 = v2-v3;
    v34 = v3-v4;
    v45 = v4-v5;
    v56 = v5-v6;
    v67 = v6-v7;
    v78 = v7-v8;
    v89 = v8-v9;
    v910 = v9-v10;
    v101 = v10-v1;
    
    % collision region normals
    v12_n = [v12(2); -v12(1)];
    v23_n = [v23(2); -v23(1)];
    v34_n = [v34(2); -v34(1)];
    v45_n = [v45(2); -v45(1)];
    v56_n = [v56(2); -v56(1)];
    v67_n = [v67(2); -v67(1)];
    v78_n = [v78(2); -v78(1)];
    v89_n = [v89(2); -v89(1)];
    v910_n = [v910(2); -v910(1)];
    v101_n = [v101(2); -v101(1)];
     
    const = [(va-v1)'*v12_n; 
            (va-v2)'*v23_n; 
            (va-v3)'*v34_n; 
            (va-v4)'*v45_n; 
            (va-v5)'*v56_n; 
            (va-v6)'*v67_n; 
            (va-v7)'*v78_n; 
            (va-v8)'*v89_n;
            (va-v9)'*v910_n;
            (va-v10)'*v101_n];     
end

%% Terminal ineq Constraints

function const = assign_terminal_ineq_constraints()
    const = assign_running_ineq_constraints();
end

function ownership = assign_terminal_ineq_ownerships()
    ownership = assign_running_ineq_ownerships();
end

%% Terminal constraints 

function const = terminal_constraint_a1(x)
    xa1 = x(1:2);
    const = left_lane_constraint(xa1);
end

function const = terminal_constraint_b1(x)
    xb1 = x(5:6);
    const = [left_lane_constraint(xb1)];
end


function const = terminal_constraint_c1(x)
    xc1 = x(9:10);
    const = right_lane_constraint(xc1);
end

function const = terminal_constraint(x)
    const = [terminal_constraint_b1(x);
             terminal_constraint_c1(x)];
end

function val = terminal_cost(~)
    val = 0;
end

function val = terminal_cost_a1(xu)
    global lane_pen lat_speed_pen;
    xa1 = xu(1:4);
    val = 0;
end

%% Running costs
function val = running_cost_a1(xu)
    global R;
    global gv_a1;
    global speed_pen lat_speed_pen lane_pen;
    xa1 = xu(1:4);
    ua1 = xu(13:14);
    val = ua1'*R*ua1 + speed_pen*speed_cost(xa1,gv_a1)+lat_speed_pen*lat_speed_cost(xa1);
end

function val = running_cost_b1(xu)
    global R;
    global gv_b1;
    global speed_pen lat_speed_pen;
    xb1 = xu(5:8);
    ub1 = xu(15:16);
    val = ub1'*R*ub1 + speed_pen*speed_cost(xb1,gv_b1)+lat_speed_pen*lat_speed_cost(xb1);
end

function val = running_cost_c1(xu)
    global R;
    global gv_c1 speed_pen lat_speed_pen;
    xc1 = xu(9:12);
    uc1 = xu(17:18);
    val = uc1'*R*uc1 + speed_pen*speed_cost(xc1,gv_c1)+lat_speed_pen*lat_speed_cost(xc1);
end

function val = running_cost_d1(xu)
    global R;
    global gv_d1 speed_pen lat_speed_pen;
    xd1 = xu(13:16);
    ud1 = xu(23:24);
    val = ud1'*R*ud1 + speed_pen*speed_cost(xd1,gv_d1)+lat_speed_pen*lat_speed_cost(xd1);
end

%% Polite costs
function val = running_cost_a1_polite(xu)
   val = running_cost_a1(xu);
end

function val = running_cost_b1_polite(xu)
    val = running_cost_b1(xu);
end

function val = running_cost_b2_polite(xu)
    val = running_cost_b2(xu);
end
%% Helpers

function val = speed_cost(x,goal)
    global linear_dyn;
    if linear_dyn
        val = (x(4)-goal)^2;
    else
        val = (x(3)-goal)^2;
    end
end

function val = lat_speed_cost(x)
    global linear_dyn;
    if linear_dyn
        val = x(3)^2;
    else
        val = 0;
    end
end

function val = either_lane_cost(x)
    global border_left border_right;
    val = ((x(1)-border_left/2)*(x(1)-border_right/2))^2;
end

function val = right_lane_cost(x)
    global border_right;
    val = (x(1)-border_right/2)^2;
end

function val = left_lane_cost(x)
    global border_left;
    val = (x(1)-border_left/2)^2;
end

function val = crosswalk_cost(x)
    global crosswalk_length crosswalk_position;
    val = (x(2)-crosswalk_position-crosswalk_length/2)^2;
end

function const = left_lane_constraint(x)
    global border_left;
    const = x(1)-border_left/2;
end

function const = right_lane_constraint(x)
    global border_right;
    const = x(1)-border_right/2;
end

function const = left_lane_ineq_constraint(x)
    global border_left half_w;
    const = x(1)-half_w-border_left;
end

function const = left_lane_double_ineq_constraint(x)
    global border_left half_w;
    const = [x(1)-half_w-border_left;
             0-(x(1)+half_w)];

end

function const = right_lane_ineq_constraint(x)
    global border_right half_w;
    const = border_right-half_w-x(1);
end

function const = avoid_constraint(x1,x2)
    global avoid_r;
    delt = x1(1:2)-x2(1:2);
    const = sqrt(delt'*delt)-avoid_r;
end

function const = crosswalk_constraint(x)
    global crosswalk_length crosswalk_position ped_w;
    const = [crosswalk_position+crosswalk_length+5-(x(2)+ped_w);
             x(2)-ped_w - crosswalk_position+5];
end

function const = crosswalk_cross_left_constraint(x)
    global border_left;
    const = [border_left-1-x(1)];
end

function pt = draw_car(center,color,alpha,angle)
    global half_l half_w;
    
    angle = angle - pi/2;
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    
    
    lower_left = R*[-half_w; -half_l];
    upper_left = R*[-half_w; +half_l];
    lower_right = R*[+half_w; -half_l];
    upper_right = R*[+half_w; +half_l];
    
   
    pt = patch('XData',[center(1)+lower_left(1),center(1)+upper_left(1),center(1)+upper_right(1),center(1)+lower_right(1),center(1)+lower_left(1)],...
               'YData',[center(2)+lower_left(2),center(2)+upper_left(2),center(2)+upper_right(2),center(2)+lower_right(2),center(2)+lower_left(2)],...
               'FaceAlpha',1,'FaceColor', color, 'EdgeColor', color);
end

function pt = draw_boundary(center,length)
    global buffer_width_wide buffer_width_narrow buffer_front buffer_behind_start buffer_front_tip;
    
    angle = 0;
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    
    lower_left = R*[-buffer_width_wide; -buffer_behind_start];
    upper_left = R*[-buffer_width_narrow; +buffer_front];
    tip = R*[0; +buffer_front_tip];
    lower_right = R*[+buffer_width_narrow; -buffer_behind_start-length];
    upper_right = R*[+buffer_width_wide; +buffer_front];
    
    color = [1,1,1];
    alpha = 0;
   
    pt = patch('XData',[center(1)+lower_left(1),center(1)+upper_left(1),center(1)+tip(1),center(1)+upper_right(1),center(1)+lower_right(1),center(1)+lower_left(1)],...
               'YData',[center(2)+lower_left(2),center(2)+upper_left(2),center(2)+tip(2),center(2)+upper_right(2),center(2)+lower_right(2),center(2)+lower_left(2)],...
               'FaceAlpha',alpha,'FaceColor', color,'EdgeColor', 'k');
%      pt = patch('XData',[center(1)+lower_left(1),center(1)+upper_left(1),center(1)+upper_right(1),center(1)+lower_right(1),center(1)+lower_left(1)],...
%                'YData',[center(2)+lower_left(2),center(2)+upper_left(2),center(2)+upper_right(2),center(2)+lower_right(2),center(2)+lower_left(2)],...
%                'FaceAlpha',alpha,'FaceColor', color,'EdgeColor', 'k');
end

function pt = draw_ped(center,color)
    global ped_w;
    pt = patch('XData',[center(1)-ped_w,center(1)-ped_w,center(1)+ped_w,center(1)+ped_w,center(1)-ped_w],...
               'YData',[center(2)-ped_w,center(2)+ped_w,center(2)+ped_w,center(2)-ped_w,center(2)-ped_w],...
               'FaceAlpha',1,'FaceColor', color);
end

function draw_lanes()
    global border_left border_right;
    plot([border_left, border_left],[-200,10000],'-k');
%     plot([border_left/2, border_left/2],[-200,200],'--k');
    plot([0, 0],[-200,10000],'--k');
%     plot([border_right/2, border_right/2],[-200,200],'--k');
    plot([border_right, border_right],[-200,10000],'-k');
end

function draw_crosswalk()
    global crosswalk_position crosswalk_length border_left border_right num_stripes;
    stripe_width = (border_right-border_left) / (2*num_stripes+1);
    loc = border_left+stripe_width;
    for i = 1:num_stripes
        draw_stripe([loc,crosswalk_position],stripe_width, crosswalk_length);
        loc = loc + 2*stripe_width;
    end
end

function draw_stripe(location, width, length)
    x = location(1);
    y = location(2);
    pt = patch('XData',[x,x,x+width,x+width,x],...
               'YData',[y,y+length,y+length,y,y],...
               'FaceAlpha',1,'FaceColor', 'w');
end


function vec = unicycle_dyn(x,u)
    global linear_dyn dt;
    vec = x(3)*cos(x(4));
    vec = [vec; x(3)*sin(x(4))];
    vec = [vec; u(1)];
    vec = [vec; u(2)];
    if linear_dyn
        vec = x(3)+dt/2*u(1);
        vec = [vec; x(4)+dt/2*u(2)];
        vec = [vec; u(1)];
        vec = [vec; u(2)];
    end
end

function next = f(xu)
    global dt;
    xa1 = xu(1:4);
    xb1 = xu(5:8);
    xc1 = xu(9:12);
%     xd1 = xu(13:16);
    
    ua1 = xu(13:14);
    ub1 = xu(15:16);
    uc1 = xu(17:18);
%     ud1 = xu(23:24);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    zc1 = xc1 + dt*unicycle_dyn(xc1,uc1);
%     zd1 = xd1 + dt*unicycle_dyn(xd1,ud1);
    next = [za1;zb1;zc1];
end