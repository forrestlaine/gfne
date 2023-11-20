% static vs dynamic example
clear; clc; 
N = 2;
T = 50;
dt = 0.1;
%%
% x0 = [5;-3;2;1];
% m{1} = 1;
% m{2} = 1;
% 
% for t = 1:T
%     F{t} = [0 1 0 dt 0 0 0;
%             0 0 1 0 dt 0 0;
%             0 0 0 1 0 dt 0;
%             0 0 0 0 1 0 dt];
%     Q{t,1} = diag([0,0,0,0,0,1,0]);
%     Q{t,2} = diag([0,0,0,0,0,0,1]);
%     H{t,1} = zeros(0,6);
%     H{t,2} = zeros(0,6);
% end
% Q{T+1,1} = 100*diag([0,0,1,0,0]);
% Q{T+1,2} = 100*blkdiag(0,[1 -1; -1 1],0,0);
% H{T+1,1} = zeros(0,5);
% H{T+1,2} = zeros(0,5);
%%
x0 = [5;-3];
m{1} = 1;
m{2} = 1;

for t = 1:T
    F{t} = [0 1 0 dt 0;
            0 0 1 0 dt];
    Q{t,1} = 1*diag([0,0,0,1,0]);
    Q{t,2} = 1*diag([0,0,0,0,1]);
    H{t,1} = zeros(0,4);
    H{t,2} = zeros(0,4);
end
Q{T+1,1} = 10*diag([0,0,1]);
Q{T+1,2} = 10*blkdiag(0,[1 -1; -1 1]);
H{T+1,1} = zeros(0,3);
H{T+1,2} = zeros(0,3);


%%
            
[X,U,L,M,P,K] = solve_ec_lq_game_super_fast(F,... % time-indexed dict of dynamics
                                                 H,... % time-indexed dict of player-indexed constraints
                                                 Q,... % time-indexed dict of player-indexed costs
                                                 N,... % number of players
                                                 T,...
                                                 x0,...
                                                 false);           
% for t= 1:T+1
%     xx(t) = X{t}(1);
%     yy(t) = X{t}(2);
% end
% 

close all;
fh = figure;

writerObj = VideoWriter('staticdynamic','MPEG-4');
writerObj.FrameRate = 20;
open(writerObj);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
set(gcf,'Color', [.81,.81,.81]);
hold on;


XX = linspace(-10,10,20);
YY = linspace(-10,10,20);

xa = [5;-3];
xb = [3; 5];
xc = [-5; 3];
xd = [-3; -5];

for t = 1:T
    ind = 1;
    for i = 1:20
        for j = 1:20
            pt = [1;XX(i);YY(j)];
            XXX(ind) = XX(i);
            YYY(ind) = YY(j);
            ctrl = K{t}*pt;
            U1(ind) = ctrl(1);
            
            U2(ind) = ctrl(2);
            ind = ind+1;
        end
    end
    q1 = quiver(XXX,YYY,U1,zeros(size(U1)),'b','Linewidth',1.25);
    q2 = quiver(XXX,YYY,zeros(size(U2)),U2,'r','Linewidth',1.25);
%     axis('off');
    axis('equal');
    axis([-10,10,-10,10]);
    
    spota = viscircles(xa',0.25,'Color', 'k');
    spotb = viscircles(xb',0.25,'Color', 'k');
    spotc = viscircles(xc',0.25,'Color', 'k');
    spotd = viscircles(xd',0.25,'Color', 'k');
    
    xa = xa + 0.1 * K{t}*[1;xa];
    xb = xb + 0.1 * K{t}*[1;xb];
    xc = xc + 0.1 * K{t}*[1;xc];
    xd = xd + 0.1 * K{t}*[1;xd];
    
    set(gca,'color',[.81 .81 .81])
    frame = getframe(fh);
    writeVideo(writerObj,frame);
    
    delete(q1);
    delete(q2);
    delete(spota);
    delete(spotb);
    delete(spotc);
    delete(spotd);
    
end

close(writerObj);

%%
figure;

set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
set(gcf,'Color', [.84,.84,.84]);
% axis('off');
axis('equal');
set(gcf, 'InvertHardCopy', 'off'); 
set(gca,'color',[.84 .84 .84])
axis([-10,10,-10,10]);
hold on;
for t= 1:T+1
    if t <= T && mod(t,5) == 0
        viscircles(X{t}',0.25,'Color', 'k');
    
        quiver(X{t}(1),X{t}(2),K{t}(1,:)*[1;X{t}],0,2,'b','Linewidth',3);
        quiver(X{t}(1),X{t}(2),0,K{t}(2,:)*[1;X{t}],2,'r','Linewidth',3);
    end
end
set(gca,'color',[.84 .84 .84])
saveas(gcf,'static_dyn.png');
