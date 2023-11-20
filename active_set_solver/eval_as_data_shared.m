function [F,H,G,Q] = eval_as_data_shared(evaluators, xval,uval,lamval,muval,gamval,T,N)
    for t = 1:T
        uu = [];
        for i = 1:N
            uu = [uu; uval{t,i}];
        end
        xu = [xval{t};uu];
        H{t} = [full(evaluators.eval_constraint{t}(xu)),...
                  full(evaluators.eval_jac_constraint{t}(xu))];
        for j = 1:numel(evaluators.eval_ineq_constraint{t})
            G{t}{j} = [full(evaluators.eval_ineq_constraint{t}{j}(xu)),...
                      full(evaluators.eval_ineq_jac_constraint{t}{j}(xu))];
        end

        F{t} = [full(evaluators.eval_pred{t}(xu))-xval{t+1}, full(evaluators.eval_jac_pred{t}(xu))];

        for i = 1:N
            cost = full(evaluators.eval_cost{t,i}(xu));
            grad = full(evaluators.eval_grad_cost{t,i}(xu));
            hess = evaluators.eval_hess_cost{t,i}(xu)+...
                   evaluators.eval_hess_constraint_mult{t,i}(xu,muval{t})+...
                   evaluators.eval_hess_dyn_mult{t,i}(xu,lamval{t,i});%+...
%                    evaluators.eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t});
            hess = full(hess);
            
            Q{t,i} = [0 grad'; 
                      grad hess];
        end
    end
    H{T+1} = [full(evaluators.eval_constraint{T+1}(xval{T+1})), full(evaluators.eval_jac_constraint{T+1}(xval{T+1}))];
    G{T+1} = {};
    for j = 1:numel(evaluators.eval_ineq_constraint{T+1})
        G{T+1}{j} = [full(evaluators.eval_ineq_constraint{T+1}{j}(xval{T+1})), full(evaluators.eval_jac_ineq_constraint{T+1}{j}(xval{T+1}))];
    end
    for i = 1:N
        cost = full(evaluators.eval_cost{T+1,i}(xval{T+1}));
        grad = full(evaluators.eval_grad_cost{T+1,i}(xval{T+1}));
        hess = evaluators.eval_hess_cost{T+1,i}(xval{T+1})+...
            evaluators.eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1});%+...
%             evaluators.eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1});
        hess = full(hess);
        Q{T+1,i} = [0 grad'; 
                    grad hess];
    end  
end