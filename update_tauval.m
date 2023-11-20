function tauval = update_tauval(params,restoration,min_comp,comp,dim,tauval,fresh_policy_count,current_residual,N)
    if params.independent_updates
        for i = 1:N
            if params.basic_tauschedule
                if private_residual{i} < tauval{i}+params.tauval_tolerance*params.tauval_decrease && tauval{i}-params.tauval_tolerance*params.tauval_decrease > params.tauval_tolerance
                    tauval{i} = params.tauval_decrease*tauval{i};
                end
            else
                if ~restoration
                    eta = min_comp{i}/(comp{i}/dim{i});
                    sigma = 0.1*min(0.05*(1-eta)/eta,2)^3;
                    tauval{i} = sigma*comp{i}/dim{i};
                end
            end
        end
    else
        if ~params.wait_for_policy_to_update_tau || fresh_policy_count >= params.fresh_policy_threshold
            if params.basic_tauschedule 
                if(current_residual < tauval{1}+params.tauval_tolerance*params.tauval_decrease)
                    for i = 1:N
                        tauval{i} = params.tauval_decrease*tauval{i};
                    end
                end
            else
                max_tau = -inf;
                for i = 1:N
                    if dim{i} > 0
                        eta = min_comp{i}/(comp{i}/dim{i});
                        sigma = 0.1*min(0.05*(1-eta)/eta,2)^3;
                        tauval{i} = sigma*comp{i}/dim{i};
                        max_tau = max(max_tau, tauval{i});
                    end
                end
                for i = 1:N
                    tauval{i} = max_tau;
                end
            end
        end
    end
end