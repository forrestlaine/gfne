%% Test fast version of ec_lq_solver
close all; clc; clear;
total_t1 = 0;
total_t2 = 0;
global_iters = 3;
for global_iter = 1:global_iters
    T = 500;
    N = 2;

    dt = 0.1;

    n = 8;
    ma = 2;
    mb = 2;
    m{1} = ma;
    m{2} = mb;

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
    H{T/2,1} = [zeros(4,1) -eye(4) eye(4) zeros(4,2)];
    H{T/2,2} = [ones(4,1)  zeros(4) eye(4) zeros(4,2)];
    % H{T/2,1} = H{T/2,1}(3:4,:);
    % H{T/2,2} = H{T/2,2}(3:4,:);

    H{T+1,1} = [zeros(4,1) -eye(4) eye(4)];
    H{T+1,2} = [ones(4,1)  eye(4) zeros(4)];
    % H{T+1,2} = zeros(0,1+n);

    % H{T+1,1} = zeros(0,1+n);
    % H{T+1,2} = zeros(0,1+n);
    avg_d = 0;
    avg_d1 = 0;

    xa(1:4,1) = [-3;-3;0;-2];
    xb(1:4,1) = [3;3;2;0];
    x0 = [xa(1:4,1);xb(1:4,1)];

    tic
    [dXa,dUa,dL,dMa,dP] = solve_ec_lq_game_super_fast(F,H,Q,N,T,x0);
    t1 = toc
    tic
    [K,k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
    thalf = toc
    [dX,dU,dL,dM,dP] = solve_ec_lq_game_d(F,H,Q,N,T,K,x0);
    t2 = toc
    total_t1 = total_t1 + t1;
    total_t2 = total_t2 + t2;
end

total_t1 = total_t1 / global_iters;
total_t2 = total_t2 / global_iters;
total_t2/total_t1

disp('what happened?');



