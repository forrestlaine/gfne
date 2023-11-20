function [xvals,solve_times,solve_iters] = mpgp_fb_game_solver_active_set_shared(f,... % dynamics function
                                                                h,... % time-indexed dict of constraint functions
                                                                hi,... % time-indexed dict of constraint ownership
                                                                g,... % time-indexed dict of inequality constraint functions
                                                                gi,... % time-indexed dict of inequality constraint ownership
                                                                l,... % time-indexed dict of player-indexed cost functions
                                                                n,... % dimension of shared state
                                                                m,... % player-indexed dict of player input dimensions
                                                                N,... % number of players
                                                                T,... % trajectory horizon
                                                                K,... % mpgp iterations
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
            working_set{t}{j}{2} = 0;
            working_set{t}{j}{3} = [];
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
        working_set{T+1}{j}{2} = 0;
        working_set{T+1}{j}{3} = [];
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

params.plot_fn = plot_fn;
% if params.debug_plot
%     plot_fn(x);
% end
 
for iter = 1:K
    
    tic;
    [dx,du,~,~,~,~,working_set_out,iters,flag] = active_set_lq_game_solver_best(F,H,hi,G,gi,Q,N,T,m,x,u,working_set,params);
    if ~flag
        for t = 1:T
            for j = 1:numel(G{t})
                working_set{t}{j}{1} = 0;
                working_set{t}{j}{2} = 0;
                working_set{t}{j}{3} = [];
            end
            allu = [];
            for i = 1:N
                u{t,i} = params.rollout_controls{i};
                allu = [allu; u{t,i}]; 
            end
            x{t+1} = full(evaluators.eval_pred{t}([x{t};allu]));                                                           
        end
        [dx,du,~,~,~,~,working_set_out,iters,flag] = active_set_lq_game_solver_best(F,H,hi,G,gi,Q,N,T,m,x,u,working_set,params);
        if ~flag
            disp('breaking mpgp early because something went wrong');
            break;
        end
    end
    
    solve_times{iter} = toc;
    solve_iters{iter} = iters;
    xvals{iter} = dx;

    if false
        p2_wrap = 0;
        p3_wrap = 0;
        p4_wrap = 0;
        dist = params.wrap_val;
        if dx{2}(6) > dx{2}(2) + dist
            p2_wrap = -1;
        elseif dx{2}(6) < dx{2}(2) - dist
            p2_wrap = 1;
        end
        if dx{2}(10) > dx{2}(2) + dist
            p3_wrap = -1;
        elseif dx{2}(10) < dx{2}(2) - dist
            p3_wrap = 1;
        end
        if dx{2}(14) > dx{2}(2) + dist
            p4_wrap = -1;
        elseif dx{2}(14) < dx{2}(2) - dist
            p4_wrap = 1;
        end


        for t = 1:T-1
            x{t} = dx{t+1};
            x{t}(6) = x{t}(6) + p2_wrap*dist*2;
            x{t}(10) = x{t}(10) + p3_wrap*dist*2;
            x{t}(14) = x{t}(14) + p4_wrap*dist*2;
            for i = 1:N
                u{t,i} = du{t+1,i};
            end
            working_set{t} = working_set_out{t+1};
        end
        allu = [];
        for i = 1:N
            u{T,i} = zeros(m{i},1);
            allu = [allu; u{T,i}];
        end
        x{T} = dx{T+1};
        x{T}(6) = x{T}(6) + p2_wrap*dist*2;
        x{T}(10) = x{T}(10) + p3_wrap*dist*2;
        x{T}(14) = x{T}(14) + p4_wrap*dist*2;
        working_set{T} = working_set_out{T+1};

        x{T+1} = F{T}*[1;x{T};allu];

        disp(['mpc iter: ' num2str(iter)]);
        if iter ==2
            disp('check');
        end

        for j = 1:numel(G{T+1})
            working_set{T+1}{j}{1} = 0;
            working_set{T+1}{j}{2} = 0;
            working_set{T+1}{j}{3} = [];
        end
    end
    
end
          
end

