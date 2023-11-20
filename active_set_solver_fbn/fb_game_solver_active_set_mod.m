function [x,working_set,total_duration] = fb_game_solver_active_set_mod(f,... % dynamics function
                                                                h,... % time-indexed dict of player-indexed constraint functions
                                                                g,... % time-indexed dict of player-indexed inequality constraint functions
                                                                gr,... % ineq regions
                                                                l,... % time-indexed dict of player-indexed cost functions
                                                                n,... % dimension of shared state
                                                                m,... % player-indexed dict of player input dimensions
                                                                N,... % number of players
                                                                T,... % time steps
                                                                x0,...% initial state
                                                                x0_rollout,...
                                                                params,...
                                                                plot_fn,...
                                                                working_set,...
                                                                initial)
all_m = 0;
for i = 1:N
    start_m{i} = all_m;
    all_m = all_m + m{i};
end
% setup
evaluators = generate_as_evaluators_mod(f,h,g,gr,l,n,m,N,T);

% Initialize solution
x{1} = x0_rollout;
for t = 1:T
    allu = [];
    for i = 1:N
        if i == 1
            u{t,i} = [15];
        else
            u{t,i} = zeros(m{i},1);
        end
        allu = [allu; u{t,i}];
        lambda{t,i} = zeros(n,1);
        mu{t,i} = zeros(size(evaluators.eval_constraint{t,i}([zeros(n+all_m,1)])));
        gam{t,i} = zeros(size(evaluators.eval_ineq_constraint{t,i}([zeros(n+all_m,1)])));
        psi{t,i} = zeros(size(all_m-m{i},1),1);
        if ~params.use_ws
            working_set{t,i} = zeros(size(gam{t,i}));
        end
    end
    x{t+1} = full(evaluators.eval_pred{t}([x{t};allu]));                                                           
end
for i = 1:N
    mu{T+1,i} = zeros(size(evaluators.eval_constraint{T+1,i}(zeros(n,1))));
    gam{T+1,i} = zeros(size(evaluators.eval_ineq_constraint{T+1,i}(zeros(n,1))));
    if ~params.use_ws
        working_set{T+1,i} = zeros(size(gam{T+1,i}));
    end
end
x{1} = x0;
if params.use_initialization
    for t = 1:T
        x{t+1} = initial{t+1};
    end
end


    
% plot_fn(x);
duration = 0;
data_duration = 0;
solve_duration = 0;
solve_iters = 0;
[F,H,G,GR,Q] = eval_as_data_mod(evaluators,x,u,lambda,mu,gam,T,N);

% if params.debug_plot
%     plot_fn(x);
% end
%     
for iter = 1:100
    if iter == 1
        unlimited = false;
    else
        unlimited = false;
    end
    tic;
    if iter > 1
        if params.linear
            [F,H,G,Q] = eval_as_data_linear_mod(evaluators,x,u,lambda,mu,gam,T,N,F,H,Q);
        else
            [F,H,G,GR,Q] = eval_as_data_mod(evaluators, x,u,lambda,mu,gam,T,N);
        end
    end
    data_time = toc;
    
    
    
    [dx,du,lambda,mu,gam,psi,working_set,iters] = active_set_lq_game_solver_mod(F,H,G,GR,Q,N,T,m,zeros(n,1),working_set,unlimited,params);
    solve_time = toc-data_time;
    solve_iters = solve_iters+iters;
    alpha = 1;
    if iter > 10 
        alpha = 1;
    end
    step = 0;
    for t = 1:T
        for i = 1:N
            u{t,i} = u{t,i} + alpha * du{t,i};
        end
        x{t+1} = x{t+1} + alpha * dx{t+1};
        step = step + norm(dx{t+1});
    end
    step
    if params.debug_plot
        plot_fn(x);
    end
    solve_duration = solve_duration + solve_time;
    data_duration = data_duration + data_time;
    duration = duration+toc;
    if step < inf
        break;
    end
end
total_duration.total = duration;
total_duration.data = data_duration;
total_duration.solve = solve_duration;
total_duration.iters = solve_iters;
end

