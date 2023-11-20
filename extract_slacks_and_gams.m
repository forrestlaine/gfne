function [alpha_dual, alpha_primal, dG, dS] = extract_slacks_and_gams(params,dX,dU,xval,uval,gamval,slackval,tauval,evaluators,T,N)

    for i = 1:N
        alpha_dual{i} = 1;
        alpha_primal{i} = 1;
    end
    for t = 1:T+1
        dxu = dX{t};
        xu = [xval{t}];
        if t < T+1
            for j = 1:N
                dxu = [dxu; dU{t,i}];
                xu = [xu; uval{t,i}];
            end
        end

        for i = 1:N
            if t < T+1
                ineq_jac = evaluators.eval_full_jac_ineq_constraint{t,i}(xu);
                ineq_val = evaluators.eval_ineq_constraint{t,i}([xval{t};uval{t,i}]);
            else
                ineq_jac = evaluators.eval_jac_ineq_constraint{T+1,i}(xval{T+1});
                ineq_val = evaluators.eval_ineq_constraint{T+1,i}(xval{T+1});
            end
            Sig = sparse(diag(gamval{t,i}./slackval{t,i}));
            dG{t,i} = full(-Sig*(ineq_jac*dxu + ineq_val - slackval{t,i} - tauval{i}*gamval{t,i}.^(-1)));
            dS{t,i} = full(Sig\(tauval{i}*slackval{t,i}.^(-1) - dG{t,i}));

            for j = 1:size(dG{t,i})
                if params.frac_to_bound*gamval{t,i}(j) + alpha_dual{i}*(dG{t,i}(j)-gamval{t,i}(j)) < 1e-12
                    alpha_dual{i} = (params.frac_to_bound*gamval{t,i}(j)-1e-12)/(gamval{t,i}(j)-dG{t,i}(j));

                end
                if params.frac_to_bound*slackval{t,i}(j) + alpha_primal{i}*dS{t,i}(j) < 1e-12
                    alpha_primal{i} = (1e-12 -params.frac_to_bound*slackval{t,i}(j))/dS{t,i}(j);
                end
            end 
        end
    end
    for i = 1:N
        if alpha_primal{i} < 0 || alpha_dual{i} < 0 || alpha_dual{i} > 1 || alpha_primal{i} > 1
            disp('VERY BAD');
        end
    end
    if ~params.independent_updates
        min_alpha_primal = inf;
        min_alpha_dual = inf;
        for i = 1:N
            min_alpha_primal = min(alpha_primal{i},min_alpha_primal);
            min_alpha_dual = min(alpha_dual{i},min_alpha_dual);
        end
        for i = 1:N
            alpha_primal{i} = min_alpha_primal;
            alpha_dual{i} = min_alpha_dual;
        end
        if params.single_alpha
            alpha_both = min(alpha_primal{1},alpha_dual{1});
            for i = 1:N
                alpha_primal{i} = alpha_both;
                alpha_dual{i} = alpha_both;
            end
        end
    end
end