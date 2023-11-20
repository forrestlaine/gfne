% Test 2d interior point method
clear; clc; close all;
import casadi.*

n = 2; m = 0; l = 2;
tau = 0.995;
eta = 0.25;
colors = ['-r','-g','-b','-k','-c','-m','-y'];

x = MX.sym('x', n); % primal variables
s = MX.sym('s', l); % slacks
y = MX.sym('y', m); % eq mults
z = MX.sym('z', l); % ineq mults

mu = MX.sym('mu', 1);
v = MX.sym('v', 1);

f = x'*diag([3,1])*x;
h = [];%[1 1]*x - 1; % [1 1]*x -1 = 0
g = [x'*x - 1; % x' * x -1 >= 0
     [1 0]*x]; % [1 0] * x >= 0
merit = f - mu * sum(log(s)) + v*norm(h, 2) + v*norm(s-g,2);

L = f;
if (m>0) 
    L = L - h'*y;
end
if (l>0)
    L = L - g'*z;
end

grad_L = gradient(L, x);
comp = z.*s - mu;
grad_merit = gradient(merit, [x;s]);

FONC = [grad_L; comp; -h; s-g];

H = jacobian(FONC, [x; s; y; z]);
% spy(H);

evalFONC = Function('grad', {x; s; y; z; mu}, {FONC});
evalH = Function('grad', {x; s; y; z; mu}, {H});
evalf = Function('grad', {x}, {f});
evalg = Function('grad', {x}, {g});
evalh = Function('grad', {x}, {h});
evalGradMerit = Function('grad', {x;s;mu;v}, {grad_merit});
evalMerit = Function('grad', {x;s;mu;v}, {merit});

xx = [-2;5];
ss = [2; 2];
yy = 0;
zz = [0; 0];

mumu = 1;
% g0 = full(evalg(xx));
% if (any(g0 < 0))
%     DISP('error, not an interior point');
% end

T = 50;
Xh = zeros(T+1,2);
Yh = zeros(T+1,m);
Zh = zeros(T+1,l);
Sh = zeros(T+1,l);
Xh(1,:) = xx;
Yh(1,:) = yy;
Zh(1,:) = zz;
Sh(1,:) = ss;

%%
NN = 50;
figure; hold on;
xxx = linspace(-5,5,NN);
yyy = linspace(-10,10,NN);

for i = 1:NN
    for j = 1:NN
        ineq = full(evalg([xxx(i),yyy(j)]));
        eq = full(evalh([xxx(i),yyy(j)]));
        for k = 1:l
            ggg{k}(i,j) = ineq(k);
        end
        for k = 1:m
            hhh{k}(i,j) = eq(k);
        end
        fff(i,j) = full(evalf([xxx(i),yyy(j)]));
    end
end

contour(xxx,yyy,fff, [0.5:0.5:3]);
labels{1} = 'f';
for k = 1:l
    contour(xxx,yyy,ggg{k}',[0,0],colors(k));
    labels{k+1} = strcat('g',num2str(k));
end
ind = l+2;
for k = 1:m
    contour(xxx,yyy,hhh{k}',[0,0],colors(ind));
    labels{ind} = strcat('h',num2str(k));
    ind = ind+1;
end

%%
tf = T;
for t = 1:T
    if (t>1)
        delete(sc);
        delete(pa);
    end
    sc = scatter(Xh(1:t,1), Xh(1:t,2));
    pa = plot(Xh(1:t,1), Xh(1:t,2), '-k');
    labels{ind} = 'Iterates';
    labels{ind+1} = 'Path';
    legend(labels)
    axis('equal')
    axis([-5,5,-5,5])
    
    St = full(evalFONC(xx, ss, yy, zz, mumu));
    e = norm(St,'Inf');
    if (e < 0.001)
        mumu = 0.15*mumu;
        if (mumu < 0.001)
            tf = t-1;
            break
        end
    end
    Ht = evalH(xx, ss, yy, zz, mumu);
    St(n+1:n+l) ./ ss;
    Ht(n+1:n+l,:) = diag(ss)\Ht(n+1:n+l,:);
    
    d = -Ht\St;
    dx = d(1:n);
    ds = d(n+1:n+l);
    dy = d(n+l+1:n+l+m);
    dz = d(n+l+m+1:end);
    
    alpha_s = 1;
    alpha_z = 1;
    
    s_violation = full((1-tau)*ss - (ss + ds));
    [worst_s_violation, idx] = max(s_violation);
    if (worst_s_violation > 0)
        alpha_s = -tau*ss(idx) / ds(idx);
    else
        alpha_s = 1;
    end
    
    z_violation = full((1-tau)*zz - (zz + dz));
    [worst_z_violation, idx] = max(z_violation);
    if (worst_z_violation > 0)
        alpha_z = -tau*zz(idx) / dz(idx);
    else
        alpha_z = 1;
    end
   
    vv = 1*norm([zz;yy],'Inf');
    merit_k = full(evalMerit(xx,ss,mumu,vv));
    grad_merit_k = full(evalGradMerit(xx,ss,mumu,vv));
    D_k = grad_merit_k'*[dx;ds] / norm([dx;ds]);
    if (full(D_k) > 0) 
        disp('BAD!')
    end
    
%     alpha_s = max(alpha_s, 0.1);
%     alpha_z = max(alpha_z, 0.1);
    
%     bti = 0;
%     bad_bt = false;
%     while (full(full(evalMerit(xx+alpha_s*dx, ss+alpha_s*ds, mumu,vv)) > merit_k + eta*alpha_s*D_k))
%         alpha_s = 0.5 * alpha_s;
%         bti = bti+1;
%         if (bti > 15)
%             bad_bt = true;
%             break
%         end
%     end
%     if (bad_bt)
%         tf = t-1;
%         disp('Bad linesearch');
%         break
%     end
    
    xx = full(xx + alpha_s*dx);
    ss = full(ss + alpha_s*ds);
    zz = full(zz + alpha_z*dz);
    yy = full(yy + alpha_z*dy);
    
    Xh(t+1,:) = xx;
    Sh(t+1,:) = ss;
    Zh(t+1,:) = zz;
    Yh(t+1,:) = yy;
end

%%

% scatter(Xh(1:tf,1), Xh(1:tf,2));
% plot(Xh(1:tf,1), Xh(1:tf,2), '-k');


