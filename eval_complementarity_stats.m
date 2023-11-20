function [min_comp, comp, dim] = eval_complementarity_stats(N,T,gamval,slackval)
    for i = 1:N
        min_comp{i} = 0;
        comp{i} = 0;
        dim{i} = 0;
    end
    for t = 1:T
        for i = 1:N
            if size(gamval{t,i},1) > 0
                comp{i} = comp{i} + gamval{t,i}'*slackval{t,i};
                min_comp{i} = min(min_comp{i}, min(gamval{t,i}.*slackval{t,i}));
            end
            dim{i} = dim{i} + size(gamval{t,i},1);
        end
    end
    for i = 1:N
        if size(gamval{T+1,i},1) > 0
            comp{i} = comp{i} + gamval{T+1,i}'*slackval{T+1,i};
            min_comp{i} = min(min_comp{i}, min(gamval{T+1,i}.*slackval{T+1,i}));
        end
        dim{i} = dim{i} + size(gamval{T+1,i},1);
    end

end