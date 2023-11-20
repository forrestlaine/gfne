%% Feedback Nash for LQ games: Riccati vs KKT
clear; clc; close all;
% Setup problem parameters

T = 2;

% dim(xa) = dim(xa) = n = 4
% dim(ua) = dim(ub) = m = 2
% loss_a = || xb || + || ua ||     = [xa;xb;ua;ub]'*Q_a*[xa;xb;ua;ub]
% loss_b = || xa-xb || + || ub ||  = [xa;xb;ua;ub]'*Q_b*[xa;xb;ua;ub]

x0 = [3; -3; 0; 0; 3; 3; 0; 0];

n = 4;
m = 2;

pc = diag([100, 100, 0, 0]);  

Q_a = [zeros(n), zeros(n), zeros(n,m), zeros(n,m);
       zeros(n), pc,   zeros(n,m), zeros(n,m);
       zeros(m,n), zeros(m,n), 1*eye(m), zeros(m);
       zeros(m,n), zeros(m,n), zeros(m), zeros(m)];

Q_b = [pc, -pc, zeros(n,m), zeros(n,m);
       -pc, pc,   zeros(n,m), zeros(n,m);
       zeros(m,n), zeros(m,n), zeros(m), zeros(m);
       zeros(m,n), zeros(m,n), zeros(m), eye(m)];

P_aT = [zeros(n), zeros(n);
       zeros(n), pc];
P_bT = [pc, -pc;
       -pc, pc];

dt = 0.1;   
a = [1 0 dt 0;
      0 1 0  dt;
      0 0 1  0;
      0 0 0  1];
b = [dt*dt/2 0;
     0 dt*dt/2;
     dt 0;
     0 dt];
               
A = [a, zeros(n); 
    zeros(n), a];
Ba = [b; zeros(n,m)];
Bb = [zeros(n,m); b];
F = [A Ba Bb];


P_a = P_aT;
P_b = P_bT;
%% Riccati version:

for t = T:-1:1

M_a = Q_a + F'*P_a*F;
M_b = Q_b + F'*P_b*F;

opt = [M_a(9:10,:); M_b(11:12,:)];
K{t} = -opt(:,9:12)\opt(:,1:8);

FK = [eye(8);
      K{t}];
P_a = FK'*M_a*FK;
P_b = FK'*M_b*FK;

end

% Now roll out
xr{1} = x0;
for t = 1:T
    ur{t} = K{t}*xr{t};
    xr{t+1} = A*xr{t} + [Ba Bb]*ur{t};
end

%% KKT system version, using K matrices from Riccati

            % states and controls    % dynamics mults    % policy mults
num_vars = (2*n+2*m)*T + (2*n)  +    (T+1)*(2*2*n)           + T*(2*m);
num_primals = (2*n+2*m)*T + (2*n);
H = zeros(num_vars);
h = zeros(num_vars,1);



ind = 0;
for t = 1:T
    ind_x = (t-1)*24;
    ind_y = (t-1)*12;
    H(ind_x+1:ind_x+12,ind_y+1:ind_y+12) = Q_a;
    H(ind_x+13:ind_x+24, ind_y+1:ind_y+12) = Q_b;
    
    ind_z = num_primals+(t-1)*20;
    
    H(ind_x+1:ind_x+12,ind_z+1:ind_z+8) = [-eye(8); zeros(4,8)];
    H(ind_x+1:ind_x+12,ind_z+17:ind_z+18) = [K{t}(3:4,:)'; zeros(2); -eye(2)];
    H(ind_x+1:ind_x+12,ind_z+21:ind_z+28) = [A'; Ba'; Bb'];
    
    H(ind_x+13:ind_x+24,ind_z+9:ind_z+16) = [-eye(8); zeros(4,8)];
    H(ind_x+13:ind_x+24,ind_z+19:ind_z+20) = [K{t}(1:2,:)'; -eye(2); zeros(2);];
    H(ind_x+13:ind_x+24,ind_z+29:ind_z+36) = [A'; Ba'; Bb'];
end

ind_x = T*24;
ind_y = T*12;
ind_z = num_primals+T*20;
H(ind_x+1:ind_x+8,ind_y+1:ind_y+8) = P_aT;
H(ind_x+9:ind_x+16,ind_y+1:ind_y+8) = P_bT;
H(ind_x+1:ind_x+8,ind_z+1:ind_z+8) = -eye(8);
H(ind_x+9:ind_x+16,ind_z+9:ind_z+16) = -eye(8);

H(ind_x+17:ind_x+24,1:8) = -eye(8);
h(ind_x+17:ind_x+24) = -x0;
for t = 1:T
    ind_xx = ind_x+24+(t-1)*8;
    ind_y = (t-1)*12;
    H(ind_xx+1:ind_xx+8,ind_y+1:ind_y+8) = A;
    H(ind_xx+1:ind_xx+8,ind_y+9:ind_y+10) = Ba;
    H(ind_xx+1:ind_xx+8,ind_y+11:ind_y+12) = Bb;
    H(ind_xx+1:ind_xx+8,ind_y+13:ind_y+20) = -eye(8);
end

sol = H\h;
for t = 1:T
    ind = (t-1)*12;
    xk{t} = sol(ind+1:ind+8);
    xd{t} = xr{t}-xk{t};
    uk{t} = sol(ind+9:ind+12);
    ud{t} = ur{t}-uk{t};
end
ind = T*12;
xk{T+1} = sol(ind+1:ind+8);
xd{T+1} = xr{T+1}-xk{T+1};

%%
if T == 2
    G = [H(2*n+1:2*n+m,:); H(2*n+2*m+2*n+m+1:4*(2*n+2*m)+4*n,:); H(4*(2*n+2*m)+4*n+2*n+1:end,:)];
    g = [h(2*n+1:2*n+m,:); h(2*n+2*m+2*n+m+1:4*(2*n+2*m)+4*n,:); h(4*(2*n+2*m)+4*n+2*n+1:end,:)];
    gg = zeros(size(g,1),2*n);
    gg(end-4*n+1:end-2*n,:) = -A;
    G = [G(:, 2*n+1:6*n+4*m) G(:, 6*n+4*m+4*n+2*m+1:end)];
    
    fb1 = G\gg;
    G2 = G;
    G2(2*m+1:2*m+2*n,2*m+1:2*m+2*n) = G(2*m+1:2*m+2*n,2*m+1:2*m+2*n) + randn(8,8);
    G2(2*m+2*n+2*m+1:2*m+2*n+2*m+2*n,2*m+1:2*m+2*n) = G(2*m+2*n+2*m+1:2*m+2*n+2*m+2*n,2*m+1:2*m+2*n) + randn(8,8);
    
    fb2 = G2\gg;
end


    