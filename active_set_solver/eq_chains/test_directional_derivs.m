% % test script
% clear; clc; close all;
% 
% n = 5;
% m = 3;
% p = 3;
% pa = 1;
% pb = 2;
% pc = p-pa-pb;
% 
% 
% Q = randn(n,n);
% Q = Q'*Q + eye(n);
% R = randn(n,m);
% 
% A = randn(pa,n);
% B = randn(pb,n);
% C = randn(pc,n);
% 
% q = randn(n,1);
% a = randn(pa,1);
% b = randn(pb,1);
% 
% Ay = randn(pa,m);
% By = randn(pb,m);
% 
% M1 = [Q -A' -B' -C'; 
%      A  zeros(pa,p);
%      B  zeros(pb,p);
%      zeros(pc,n+pa+pb) eye(pc)];
% N1 = [R; Ay; By; zeros(pc,m)];
% m1 = [q; a; b; zeros(pc,1)];
% 
% M2 = [Q -A' -B' -C';
%      A zeros(pa,p);
%      zeros(pb+pc, n+pa) eye(pb+pc)];
% N2 = [R; Ay; zeros(pb+pc,m)];
% m2 = [q; a; zeros(pb+pc,1)];
% 
% K1 = M1\[N1 m1];
% K2 = M2\[N2 m2];
% 
% Kd = K1-K2;
% 
% AAA = Kd(:,1:end-1);
% bbb = Kd(:,end);
% 
% % y = -(AAA'*AAA)\AAA'*bbb;
% y = -pinv(AAA)*bbb;
% 
% 
% B*(K2(1:n,:)*[y;1])-By*y-b;
% K1(n+pa+1:n+pa+pb,:)*[y;1];
% 
% grads = [ B*K2(1:n,:) - [By b];
%           -K1(n+pa+1:n+pa+pb,:)];
% 
% ry = randn(m,1);
% eps = 1e-5;
% yy = y+eps*ry;
% % grads*[y;1]
% % grads*[yy;1]
% rank(grads)

%%
% clear; clc;
% 
% Q = randn(3,3);
% Q = Q'*Q + eye(3);
% 
% l = 2;
% A = randn(l,3);
% A2 = [A;3*A];
% 
% M = [Q -A';
%      A zeros(l)];
% Nx = randn(3,2);
% Ny = randn(l,2);
% N = [Nx; Ny];
% 
% M2 = [Q -A2';
%      A2 zeros(2*l)];
% N2 = [N; 3*Ny];
% 
% q = randn(3,1);
% b = randn(l,1);
% b2 = [b;b];
% 
% m = [q; b];
% m2 = [q; b2];
% 
% sol = M\N
% sol2 = pinv(M2)*N2

%%
% 
% A = [-.5 -1; -.5 1];
% M = [eye(2) -A'; A zeros(2)];
% 
% Mi = inv(M);
% Si = Mi(3:4,3:4);
% 
% 
% N = 50;
% T =  linspace(-pi,pi,N);
% C = zeros(2,N);
% L = zeros(2,N);
% Z = zeros(2,N);
% for t = 1:N
%     z = [cos(T(t)); sin(T(t))];
%     Z(:,t) = z;
%     C(:,t) = A*z;
%     L(:,t) = Si*C(:,t);
% end
% 
% c1_active = C(1,:) < 0;
% c2_active = C(2,:) < 0;
% l1_active = L(1,:) < 0;
% l2_active = L(2,:) < 0;
% figure;
% hold on;
% plot(T, c1_active,'linewidth',5);
% plot(T, c2_active,'linewidth',5);
% plot(T, l1_active,'linewidth',5);
% plot(T, l2_active,'linewidth',5);
% plot([-2.6779,-2.6679],[-1,2],':k');
% plot([2.6779,2.6679],[-1,2],':k');
% plot([-1.1071,-1.1071],[-1,2],':k');
% plot([1.1071,1.1071],[-1,2],':k');
% legend('c1', 'c2', 'l1', 'l2');

% %% 
% clear; clc; close all;
% 
% A = [-3 -sqrt(3) 1;
%      3 -sqrt(3) 1;
%      0 2*sqrt(3) 1];
% 
% M = [eye(3) -A';
%     A zeros(3)];
% 
% Mi = inv(M);
% Si = Mi(4:6,4:6);
% 
% N = 50000;
% 
% for a = 1:2
%     for b = 1:2
%         for c = 1:2
%             Z{a,b,c} = zeros(3,0);
%             C{a,b,c} = zeros(3,0);
%             S{a,b,c} = zeros(3,0);
%         end
%     end
% end
% mag = 10;
% options = optimoptions('quadprog', 'Display', 'off');
% for n = 1:N
%     z = mag*rand(3,1)-mag/2;
%     c = A*z;
%     lam = -Si*c;
%     [sol,fval,exitflag,output,lambda] = quadprog(eye(3), -z, -A, zeros(3,1),[],[],[],[],[], options);
%     lambda = lambda.ineqlin;
%     code = (lambda > 1e-6)+1;
%     
%     Z{code(1),code(2),code(3)} = [Z{code(1),code(2),code(3)} z];
%     S{code(1),code(2),code(3)} = [S{code(1),code(2),code(3)} lam];
%     C{code(1),code(2),code(3)} = [C{code(1),code(2),code(3)} c];
% end
% 
% 
% figure;
% hold on;
% count = 0;
% for a = 1:2
%     for b = 1:2
%         for c = 1:2
%             
%             XYZ = Z{a,b,c};
%             
%             color = 2-[a,b,c];
%             if all(color == [0 0 0])
%                 color = [1,0,0]
%                 scatter3(XYZ(1,:),XYZ(2,:),XYZ(3,:),150,color);
%             end
% 
%         end
%     end
% end


Q = diag([1,1,0]);
q = -[-1; 1; 0];
Aeq = [-1 1 -1];
beq = [0];
A = -[0 1 0;
     0 0 1];
b = [0;0];

[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = quadprog(Q,q,A,b,Aeq,beq,[],[],randn(3,1));
X