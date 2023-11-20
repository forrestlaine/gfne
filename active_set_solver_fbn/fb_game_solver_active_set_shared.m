function [x,working_set,total_duration] = fb_game_solver_active_set_shared(f,... % dynamics function
                                                                h,... % time-indexed dict of constraint functions
                                                                hi,... % time-indexed dict of constraint ownership
                                                                g,... % time-indexed dict of inequality constraint functions
                                                                gi,... % time-indexed dict of inequality constraint ownership
                                                                l,... % time-indexed dict of player-indexed cost functions
                                                                n,... % dimension of shared state
                                                                m,... % player-indexed dict of player input dimensions
                                                                N,... % number of players
                                                                T,... % time steps
                                                                x0,...% initial state
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
evaluators = generate_as_evaluators_shared(f,h,hi,g,gi,l,n,m,N,T);
x{1} = zeros(n,1);
for t = 1:T
    allu = [];
    for i = 1:N
        u{t,i} = zeros(m{i},1);
        allu = [allu; u{t,i}];
        lambda{t,i} = zeros(n,1);
        psi{t,i} = zeros(size(all_m-m{i},1),1);
    end
    mu{t} = zeros(size(evaluators.eval_constraint{t}([zeros(n+all_m,1)])));
    gam{t} = zeros(numel(g{t}),1);
    
    if numel(mu{t}) ~= numel(hi{t}) || numel(gam{t}) ~= numel(gi{t})
        disp('Check constraint ownerships!');
    end
        
    if ~params.use_ws
        working_set{t} = {};
        for j = 1:size(gam{t},1)
            working_set{t}{j}{1} = 0;
            working_set{t}{j}{2} = [];
        end
    end
    x{t+1} = zeros(n,1);                                                           
end
mu{T+1} = zeros(size(evaluators.eval_constraint{T+1}(zeros(n,1))));
gam{T+1} = zeros(numel(g{T+1}),1);

if numel(mu{T+1}) ~= numel(hi{T+1}) || numel(gam{T+1}) ~= numel(gi{T+1})
        disp('Check constraint ownerships!');
end
if ~params.use_ws
    working_set{T+1} = {};
    for j = 1:size(gam{T+1},1)
        working_set{T+1}{j}{1} = 0;
        working_set{T+1}{j}{2} = [];
    end
end

[F,H,G,Q] = eval_as_data_shared(evaluators,x,u,lambda,mu,gam,T,N);



% Initialize solution
x{1} = x0;
for t = 1:T
    allu = [];
    for i = 1:N
        u{t,i} = params.rollout_controls{i};
        allu = [allu; u{t,i}]; 
    end
    x{t+1} = full(evaluators.eval_pred{t}([x{t};allu]));                                                           
end



    
% plot_fn(x);
duration = 0;
data_duration = 0;
solve_duration = 0;
solve_iters = 0;

params.plot_fn = plot_fn;
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
        if false
            [F,H,G,Q] = eval_as_data_linear_mod(evaluators,x,u,lambda,mu,gam,T,N,F,H,Q);
        else
            [F,H,G,GR,Q] = eval_as_data_shared(evaluators, x,u,lambda,mu,gam,T,N);
        end
    end
    data_time = toc;
    
    
    
    [dx,du,lambda,mu,gam,psi,working_set,iters] = active_set_lq_game_solver_best(F,H,hi,G,gi,Q,N,T,m,x,u,working_set,params);
    solve_time = toc-data_time;
    solve_iters = solve_iters+iters;
    alpha = 1;
    if iter > 10 
        alpha = 1;
    end
    step = 0;
    for t = 1:T
        for i = 1:N
%             u{t,i} = u{t,i} + alpha * du{t,i};
            u{t,i} = du{t,i};
        end
%         x{t+1} = x{t+1} + alpha * dx{t+1};
        x{t+1} = dx{t+1};
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

