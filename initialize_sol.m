function [xval,...
          uval,...
          lamval,...
          muval,...
          slackval,...
          gamval,...
          devval,...
          devopt,...
          psival,...
          tauval] = initialize_sol(params,z0,T,N,n,m,all_m,evaluators,initialization)
      
    for i = 1:N
        tauval{i} = params.tauval_init{i};
    end

    if params.use_initialization % use provided initialization
        xval{1} = z0;
        for t = 1:T-params.init_step
            for i = 1:N
                uval{t,i} = initialization{2}{t,i};
                lamval{t,i} = initialization{3}{t,i};
                muval{t,i} = initialization{4}{t,i};
                slackval{t,i} = max(full(evaluators.eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),tauval{i});
                gamval{t,i} = tauval{i}./slackval{t,i};
                devval{t,i} = zeros(m{i},1);
                devopt{t,i} = zeros(m{i},1);
                psival{t,i} = initialization{5}{t,i};
                
            end
            xval{t+1} = initialization{1}{t+1};
        end
        for t = T-params.init_step+1:T
            xu = xval{t};
            for i = 1:N
                uval{t,i} = zeros(m{i},1);
                devval{t,i} = zeros(m{i},1);
                devopt{t,i} = zeros(m{i},1);
                lamval{t,i} = zeros(n,1);
                muval{t,i} = zeros(size(evaluators.eval_constraint{t,i}([xval{t};uval{t,i}]),1),1);
                slackval{t,i} = max(full(evaluators.eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),tauval{i});
                gamval{t,i} = tauval{i}./slackval{t,i};
                
                psival{t,i} = zeros(all_m,1);
                xu = [xu; uval{t,i}];
            end
            xval{t+1} = full(evaluators.eval_pred{t}(xu));
        end
        for i = 1:N
            muval{T+1,i} = zeros(size(evaluators.eval_constraint{T+1,i}(xval{T+1}),1),1);
            slackval{T+1,i} = max(full(evaluators.eval_ineq_constraint{T+1,i}(xval{T+1})),tauval{i});
            gamval{T+1,i} = tauval{i}./slackval{T+1,i};
        end
                
    else % initialize from zero control and zero multipliers
        if params.use_rollout_x0
            xval{1} = params.rollout_x0;
        else
            xval{1} = z0;
        end
        for t = 1:T
            xu = [xval{t}];
            for i = 1:N
                uval{t,i} = zeros(m{i},1);
                devval{t,i} = zeros(m{i},1);
                devopt{t,i} = zeros(m{i},1);
                lamval{t,i} = zeros(n,1);
                muval{t,i} = zeros(size(evaluators.eval_constraint{t,i}([xval{t};uval{t,i}]),1),1);
                slackval{t,i} = max(full(evaluators.eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),tauval{i});
                gamval{t,i} = tauval{i}./slackval{t,i};
                psival{t,i} = zeros(all_m,1);
                xu = [xu; uval{t,i}];
            end
            xval{t+1} = full(evaluators.eval_pred{t}(xu));
        end
        xval{1} = z0;
    end
    
    for i = 1:N
        muval{T+1,i} = zeros(size(evaluators.eval_constraint{T+1,i}(xval{T+1}),1),1);
        slackval{T+1,i} = max(full(evaluators.eval_ineq_constraint{T+1,i}(xval{T+1})),tauval{i});
        gamval{T+1,i} = tauval{i}./slackval{T+1,i};
    end
end