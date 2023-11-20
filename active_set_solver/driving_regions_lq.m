%% Is there a car?
clear; clc; close all;

global dt;
dt = 0.1;

global linear_dyn;
linear_dyn = true;

global border_left;
global border_right;
global lane_width;


lane_width = 4;
border_left = -lane_width;
border_right = lane_width;

global center_lane_pen speed_pen lane_pen lat_speed_pen;
center_lane_pen = 0;

lane_pen = 20000;
global R;
R = diag([1,1]);

global car_r ped_r avoid_r half_l half_w;
global buffer_length;
global buffer_width;
half_l = 1.5;
half_w = 0.75;
car_r = sqrt(half_l^2+half_w^2);
ped_r = 2;
avoid_r = (car_r+car_r);
buffer_length = 10;
buffer_width = 2.5;

global ped_w;
ped_w = 0.5;

global crosswalk_position crosswalk_length num_stripes;
crosswalk_position = 0;
crosswalk_length = 2;
num_stripes = 6;

global polite_pen;
polite_pen = 5;
speed_pen = .5;
lat_speed_pen = 100;

global gv_a1 gv_b1 gv_c1;



global T
T = 100;
needed_speed = -10/(T*dt);

gv_a1 = 25;
gv_b1 = 24;
gv_c1 = 22;

x0 = [border_right/2; 0;0; gv_a1;
      border_left/2; -10;0;gv_b1;
      border_right/2; 30;0;gv_c1];

rollout_controls{1} = [0; -.2];
rollout_controls{2} = [0; 0];
rollout_controls{3} = [0; 0];
  
m{1} = 2;
m{2} = 2;
m{3} = 2;

full_scope{1} = true;
full_scope{2} = true;
full_scope{3} = true;

N = 3;
n = N*4;

for t = 1:T
    
    lf{t,1} = @running_cost_a1;
    lf{t,2} = @running_cost_b1;
    lf{t,3} = @running_cost_c1;
  
    gf{t} = assign_running_ineq_constraints();
    gi{t} = assign_running_ineq_ownerships();
    
    h{t} = @empty_constraint;
    hi{t} = {};

end
lf{T+1,1} = @terminal_cost_a1;
lf{T+1,2} = @empty_cost;
lf{T+1,3} = @empty_cost;

gf{T+1} = assign_terminal_ineq_constraints();
gi{T+1} = assign_terminal_ineq_ownerships();

h{T+1} = @terminal_constraint;
hi{T+1}{1} = 2;
hi{T+1}{2} = 3;
% hi{T+1}{3} = 3;


%% Run 

params.use_ws = false;
params.open_loop = false;
params.use_initialization = false;
params.debug_plot = true;
params.linear = false;
params.full_scope = full_scope;
params.print_frequency = 10;
params.rollout_controls = rollout_controls;

[xval,ws,durs] = fb_game_solver_active_set_shared(@f, h, hi, gf, gi, lf, n, m, N, T, x0, params, @plot_solution);
solution_set{1} = xval;
durations{1} = durs;

%%
plot_solution(xval,0);
 
%% Plotting

function plot_solution(xval,x0)
    make_vid = false;
    make_frames = false;
    close all;
    fh = figure('Position', [300 300 150 750]); 
    if make_vid
        writerObj = VideoWriter('gfne_paper_polite','MPEG-4');
        writerObj.FrameRate = 20;
        open(writerObj);
%         axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        set(gcf,'Color', [.81,.81,.81]);
    end
    hold on;
    global car_r T;
    
    draw_lanes();

    for t = 1:T+1
        xx(t,:) = xval{t};
        
%         spot_b1 = viscircles(xx(t,5:6),car_r,'Color', [1,.7,0],'LineStyle',':');
        border_b1 = draw_boundary(xx(t,5:6));
        car_b1 = draw_car(xx(t,5:6),[1,.7,0],1,pi/2);
        
                
%         spot_a1 = viscircles(xx(t,1:2),car_r,'Color', 'b','LineStyle',':');
        car_a1 = draw_car(xx(t,1:2),'b',1,pi/2);
        
%         spot_c1 = viscircles(xx(t,9:10),car_r,'Color', 'g','LineStyle',':');
        border_c1 = draw_boundary(xx(t,9:10));
        car_c1 = draw_car(xx(t,9:10),'g',1,pi/2);
        
        
        axis([-4,+4,xx(t,10)-65,xx(t,10)+20]);
        axis('equal');
        if make_vid
            frame = getframe(fh);
            writeVideo(writerObj,frame);
        elseif make_frames
            axis([-4,+4,xx(t,10)-65,xx(t,10)+20]);
            axis('equal');
            if mod(t,floor(T/5)) == 1
                name = ['gfne_polite_frame_' num2str(t) '.png'];
                saveas(gcf,name)
            end
        else
            pause(.01);
        end
%         delete(spot_a1);
        delete(car_a1);
%         delete(spot_b1);
        delete(car_b1);
%         delete(spot_c1);
        delete(car_c1);
        delete(border_b1);
        delete(border_c1);
        
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
    % collision zone of xc
    global buffer_length;
    global buffer_width;
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    const = [xa1(2) - (xc1(2)+buffer_length);
             (xc1(2)-buffer_length) - xa1(2);
             (xc1(1)-buffer_width)-xa1(1)];  
end

function const = running_ineq_constraint_a2(x)
    % collision zone of xb
    global buffer_length;
    global buffer_width;
    xa1 = x(1:2);
    xb1 = x(5:6);
    xc1 = x(9:10);
    const = [xa1(2) - (xb1(2)+buffer_length);
             (xb1(2)-buffer_length) - xa1(2);
             xa1(1) - (xb1(1)+buffer_width)];  
end

function const = assign_running_ineq_constraints()
    const{1} = @running_ineq_constraint_a1;
    const{2} = @running_ineq_constraint_a2;  
end

function ownerships = assign_running_ineq_ownerships()
    ownerships{1} = 1;
    ownerships{2} = 1;
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
    val = lane_pen*left_lane_cost(xa1) + lat_speed_pen*lat_speed_cost(xa1);
end

%% Running costs
function val = running_cost_a1(xu)
    global R;
    global gv_a1;
    global speed_pen lat_speed_pen lane_pen;
    xa1 = xu(1:4);
    ua1 = xu(13:14);
    ub1 = xu(15:16);
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
               'FaceAlpha',alpha,'FaceColor', color);
end

function pt = draw_boundary(center)
    global half_l half_w buffer_width buffer_length;
    
    angle = 0;
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    
    width = buffer_width - half_w;
    length = buffer_length - half_l;
    
    lower_left = R*[-width; -length];
    upper_left = R*[-width; +length];
    lower_right = R*[+width; -length];
    upper_right = R*[+width; +length];
    
    color = [1,1,1];
    alpha = 1;
   
    pt = patch('XData',[center(1)+lower_left(1),center(1)+upper_left(1),center(1)+upper_right(1),center(1)+lower_right(1),center(1)+lower_left(1)],...
               'YData',[center(2)+lower_left(2),center(2)+upper_left(2),center(2)+upper_right(2),center(2)+lower_right(2),center(2)+lower_left(2)],...
               'FaceAlpha',alpha,'FaceColor', color);
end

function pt = draw_ped(center,color)
    global ped_w;
    pt = patch('XData',[center(1)-ped_w,center(1)-ped_w,center(1)+ped_w,center(1)+ped_w,center(1)-ped_w],...
               'YData',[center(2)-ped_w,center(2)+ped_w,center(2)+ped_w,center(2)-ped_w,center(2)-ped_w],...
               'FaceAlpha',1,'FaceColor', color);
end

function draw_lanes()
    global border_left border_right;
    plot([border_left, border_left],[-200,800],'-k');
%     plot([border_left/2, border_left/2],[-200,200],'--k');
    plot([0, 0],[-200,800],'--k');
%     plot([border_right/2, border_right/2],[-200,200],'--k');
    plot([border_right, border_right],[-200,800],'-k');
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
    
    ua1 = xu(13:14);
    ub1 = xu(15:16);
    uc1 = xu(17:18);
    
    za1 = xa1 + dt*unicycle_dyn(xa1,ua1);
    zb1 = xb1 + dt*unicycle_dyn(xb1,ub1);
    zc1 = xc1 + dt*unicycle_dyn(xc1,uc1);
    next = [za1;zb1;zc1];
end