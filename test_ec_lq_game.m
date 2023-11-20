%% Test ec_lq_solver
close all; clc; clear;
T = 100;
N = 2;

dt = 0.1;

n = 8;
ma = 2;
mb = 2;

Aa = randn(4);
Aa(1,3) = dt;
Aa(2,4) = dt;
Ba = randn(4,2);
Ba(1,1) = dt*dt/2;
Ba(3,1) = dt;
Ba(2,2) = dt*dt/2;
Ba(4,2) = dt;

Qa = zeros(1+n+ma+mb);
Qa(1+n+1:1+n+ma,1+n+1:1+n+ma) = eye(ma);
Qb = zeros(1+n+ma+mb);
Qb(1+n+ma+1:1+n+ma+mb,1+n+ma+1:1+n+ma+mb) = eye(mb);

A = blkdiag(Aa,Aa);
B = blkdiag(Ba,Ba);
c = randn(n,1);
Ft = [c A B];

H = cell(T,N);
Q = cell(T,N);

for t = 1:T
    F{t} = Ft;
    H{t,1} = zeros(0,1+n+ma);
    H{t,2} = zeros(0,1+n+mb);
    
    Q{t,1} = Qa;
    Q{t,2} = Qb;
end
Q{T+1,1} = zeros(1+n);
Q{T+1,2} = zeros(1+n);
% Q{T+1,2}(2:5,2:5) = eye(4);
H{T+1,1} = [zeros(4,1) -eye(4) eye(4)];
H{T+1,2} = [ones(4,1) eye(4) zeros(4)];
% H{T+1,2} = zeros(0,1+n);

% H{T+1,1} = zeros(0,1+n);
% H{T+1,2} = zeros(0,1+n);
avg_d = 0;
avg_d1 = 0;

xa(1:4,1) = [-3;-3;0;-2];
xb(1:4,1) = [3;3;2;0];
x0 = [xa(1:4,1);xb(1:4,1)];

num_iters = 1;
for it = 1:num_iters
    tic;
    [K,k] = solve_ec_lq_game_ol(F,H,Q,N,T);
    % d1 = toc;
    [dX,dU,dL,dM,dP] = solve_ec_lq_game_d(F,H,Q,N,T,K,x0);
    duration = toc;
    % duration = toc-d1;
%     [K2,k2] = solve_ec_lq_game(F,H,Q,N,T);
    duration2 = toc-duration;
    avg_d = avg_d + duration;
    avg_d1 = avg_d1 + duration2;
end
avg_d = avg_d / num_iters;
avg_d1 = avg_d1 / num_iters;

% [K3,k3] = solve_ec_lq_game_d(F,H,Q,N,T,K2,k2);




xa2(1:4,1) = xa(1:4,1);
xb2(1:4,1) = xb(1:4,1);

xa3(1:4,1) = xa(1:4,1);
xb3(1:4,1) = xb(1:4,1);



costa1 = 0;
costa2 = 0;
costa3 = 0;
costb1 = 0;
costb2 = 0;
costb3 = 0;

dbr = 1;
for t = 1:T
    ua = K{t,1}*[xa(:,t);xb(:,t)] + k{t,1};
    ub = K{t,2}*[xa(:,t);xb(:,t)] + k{t,2};
%     ua = [1;1];
    xa(:,t+1) = Aa*xa(:,t) + Ba*(ua + dbr*randn(ma,1));
    xb(:,t+1) = Aa*xb(:,t) + Ba*(ub + dbr*randn(ma,1));
    
    costa1 = costa1 + ua'*ua;
    costb1 = costb1 + ub'*ub;
    
    ua2 = K2{t,1}*[xa2(:,t);xb2(:,t)] + k2{t,1};
    ub2 = K2{t,2}*[xa2(:,t);xb2(:,t)] + k{t,2};
%     ua2 = [1;1];
    xa2(:,t+1) = Aa*xa2(:,t) + Ba*(ua2 + dbr*randn(ma,1));
    xb2(:,t+1) = Aa*xb2(:,t) + Ba*(ub2 + dbr*randn(ma,1));
    
    ua3 = K3{t,1}*[xa3(:,t);xb3(:,t)] + k3{t,1};% + 1*randn(ma,1);
    ub3 = K3{t,2}*[xa3(:,t);xb3(:,t)] + k3{t,2};% + 1*randn(ma,1);
    xa3(:,t+1) = Aa*xa3(:,t) + Ba*(ua3 + dbr*randn(ma,1));
    xb3(:,t+1) = Aa*xb3(:,t) + Ba*(ub3 + dbr*randn(ma,1));
        
    costa2 = costa2 + ua2'*ua2;
    costb2 = costb2 + ub2'*ub2;
%     
    costa3 = costa3 + ua3'*ua3;
    costb3 = costb3 + ub3'*ub3;
end

% ind = 0;
% for t = 1:T-1
%     xa3(:,t+1) = full_sol(ind+ma+mb+ma+mb+n+n+1:ind+ma+mb+ma+mb+n+n+4);
%     xb3(:,t+1) = full_sol(ind+ma+mb+ma+mb+n+n+5:ind+ma+mb+ma+mb+n+n+8);
%     ind = ind+ma+mb+ma+mb+n+n+n;
% end
% ind = ind-ma-mb;
% xa3(:,T+1) = full_sol(ind+ma+mb+ma+mb+n+n+1:ind+ma+mb+ma+mb+n+n+4);
% xb3(:,T+1) = full_sol(ind+ma+mb+ma+mb+n+n+5:ind+ma+mb+ma+mb+n+n+8);

% writerObj = VideoWriter('ec_lqgame_noise','MPEG-4');
% writerObj.FrameRate = 10;
% open(writerObj);
% % Z = peaks; surf(Z); 
% axis tight
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
% hold on;

figure;
hold on;
for t = 1:T+1
    spot_a = plot(xa(1,t), xa(2,t),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    spot_b = plot(xb(1,t),xb(2,t),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
      axis([-10 10 -10 10]);
%     frame = getframe;
%     writeVideo(writerObj,frame);
    pause(dt);
%     delete(spot_a);
%     delete(spot_b);
end
for t = 1:T+1
    spot_a = plot(xa3(1,t), xa3(2,t),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    spot_b = plot(xb3(1,t),xb3(2,t),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    axis([-10 10 -10 10]);
%     frame = getframe;
%     writeVideo(writerObj,frame);
    pause(dt);
%     delete(spot_a);
%     delete(spot_b);
end
for t = 1:T
    spot_a = plot(dX{t}(1), dX{t}(2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'g');
    spot_b = plot(dX{t}(5),dX{t}(6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'g');
    axis([-10 10 -10 10]);
%     frame = getframe;
%     writeVideo(writerObj,frame);
    pause(dt);
%     delete(spot_a);
%     delete(spot_b);
end

% close(writerObj);
% plot(xa(1,:),xa(2,:),'-b');
% plot(xb(1,:),xb(2,:),'-r');