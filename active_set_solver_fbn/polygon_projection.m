% Polygon projection
clear; clc; close all;
% dt = 0.1;
% T = 20;
% 
% a = [1 0 dt 0;
%      0 1  0 dt;
%      0 0  1 0;
%      0 0  0 1];
% b = [dt*dt/2 0;
%       0 dt*dt/2;
%       dt 0;
%       0 dt];
%   
% term_pen = 1000;
% bound = 3;
%   
% QR = zeros(6);
% QR(5:6,5:6) = 0.1*eye(2);
% for t = 1:T
%     QQ{t} = QR;
% end
% QQ{T+1} = term_pen*eye(4);
% QQ = blkdiag(QQ{:});
% qq = zeros(T*6+4,1);
% 
% xg = [2003;1;-4;-9];
% x0 = ones(4,1);
% qq(end-3:end) = -term_pen*xg;
%   
% 
% F = zeros((T+1)*4, (T+1)*4 + T*2);
% f = zeros((T+1)*4,1);
% 
% F(1:4,1:4) = eye(4);
% f(1:4) = x0;
% 
% i1 = 4;
% i2 = 0;
% for t = 1:T
%     F(i1+1:i1+4,i2+1:i2+4) = -a;
%     F(i1+1:i1+4,i2+5:i2+6) = -b;
%     F(i1+1:i1+4,i2+7:i2+10) = eye(4);
%     i1 = i1+4;
%     i2 = i2+4+2;
% end
% 
% A = zeros(0,T*6+4);
% b = [];
% i1 = 0;
% for t = 1:T
%     d = zeros(4,T*6+4);
%     d(:,i1+5:i1+6) = [1 1;
%                      -1 1;
%                      -1 -1;
%                       1 -1];
%     A = [A; d];
%     b = [b; ones(4,1)*bound];
%     i1 = i1+6;
% end
% 
% % solution 1   
% tic;
% X2 = quadprog(QQ,qq,A,b,F,f);
% t1 = toc
% 
% %%
% tic;
% [U,S,V] = svd(F);
% rk = rank(S);
% U1 = U(:,1:rk);
% V1 = V(:,1:rk);
% S1 = S(1:rk,1:rk);
% V2 = V(:,rk+1:end);
% 
% % x = V2*z2 + S1\U1'*f
% 
% % z2 = V2'*X;
% % XX = V2*z2 + V1*(S1\U1'*f);
% 
% M = V2;
% m = V1*(S1\U1'*f);
% 
% Qh = M'*QQ*M;
% Qh = (Qh+Qh')/2;
% qh = M'*(QQ*m + qq);
% Ah = A*M;
% bh = b-A*m;
% 
% % Z = quadprog(Qh,qh,Ah,bh);
% dim = size(Ah,2);
% 
% P = sqrtm(Qh);
% 
% Qg = eye(dim);
% qg = (P')\qh;
% Ag = Ah/P;
% bg = bh;
% 
% for l = 1:size(Ag,1)
%     a = Ag(l,:);
%     n = norm(a);
%     Ag(l,:) = Ag(l,:) / n;
%     bg(l) = bg(l) / n;
% end   
% U = quadprog(Qg,qg,Ag,bg);
% r = Ag*U-bg;
% a_inds = find(r >= -1e-6);
% 
% Z = P\U;
% X = M*Z+m;
% 
% Uh = -qg;
% rh = Ag*Uh-bg;
% a2_inds = find(rh >= 0);
% 
% both_inds = intersect(a_inds,a2_inds);
% 
% t2 = toc;
% 
% if size(a2_inds,1) == size(both_inds,1)
%     disp('all initial violations are active');
% else
%     disp('not true');
% end

A = [0 1 1;
     0 0 -1;
     1 1 1;
     -1 1 1;
     1 -1 1;
     -1 -1 1];
b = [.2;2;1;1;3;2];

for i = 1:size(A,1)
    nm = norm(A(i,:));
    A(i,:) = A(i,:)/nm;
    b(i) = b(i) / nm;
end


figure;
hold on;
[V,nr] = con2vert(A,b);

[k2,av2] = convhull(V(:,1),V(:,2),V(:,3),'Simplify',true);

trisurf(k2,V(:,1),V(:,2),V(:,3),'FaceColor','cyan')
axis equal;

% for i = 1:size(V,1)
%     vtxs = [V(i,:)];
%     for j = 1:size(A,1)
%         if abs(A(j,:)*V(i,:)' - b(j)) < 1e-5
%             vtxs = [vtxs; V(i,:)+20*A(j,:)];
% %             mArrow3(V(i,:),V(i,:)+20*A(j,:));
%         end
%     end 
%     [k2,av2] = convhull(vtxs(:,1),vtxs(:,2),vtxs(:,3),'Simplify',true);
%     trisurf(k2,vtxs(:,1),vtxs(:,2),vtxs(:,3),'FaceColor','r','FaceAlpha',0.1);
% end

% for i = 1:size(A,1)
%     vtxs = [];
%     for j = 1:size(V,1)
%         if abs(A(i,:)*V(j,:)' - b(i)) < 1e-5
%             vtxs = [vtxs; V(j,:); V(j,:)+20*A(i,:)];
%         end
%     end
%     [k2,av2] = convhull(vtxs(:,1),vtxs(:,2),vtxs(:,3),'Simplify',true);
%     trisurf(k2,vtxs(:,1),vtxs(:,2),vtxs(:,3),'FaceColor','g','FaceAlpha',0.1);
% end

for k = 1:1000
    pt = rand(3,1);
    pt = pt*20 - 10;
    
    pt2 = quadprog(eye(3),-pt,A,b);
    
    cons2 = A*pt2 - b;
    inds2 = find(cons2 > -1e-6);
    
    cons = A*pt-b;
    inds = find(cons >= 0);
    
    [m,i] = max(cons);
    if ~ismember(i,inds2) && norm(pt-pt2) > 1e-3

        scatter3(pt2(1),pt2(2),pt2(3),15,'r','filled');
%         figure;
        hold on;
        scatter3(pt(1),pt(2),pt(3),15,'g','filled');
        for j = 1:size(A,1)
            a = A(j,:);
            
            xy = [-10;-10];
            z = (b(j)-a(1:2)*xy)/a(3);
            bl = [xy' z];
            
            xy = [-10;10];
            z = (b(j)-a(1:2)*xy)/a(3);
            tl = [xy' z];
            
            xy = [10;10];
            z = (b(j)-a(1:2)*xy)/a(3);
            tr = [xy' z];
            
            xy = [10;-10];
            z = (b(j)-a(1:2)*xy)/a(3);
            br = [xy' z];
            
            corners = [bl;tl;tr;br;bl];
            if j == i
                patch(corners(:,1),corners(:,2),corners(:,3),'b','FaceAlpha',0.1);
            elseif ismember(j,inds2)
                patch(corners(:,1),corners(:,2),corners(:,3),'g','FaceAlpha',0.1);
% %             else
% %                 patch(corners(:,1),corners(:,2),corners(:,3),'r','FaceAlpha',0.1);
            end
        end
        disp('what');
    end
end

