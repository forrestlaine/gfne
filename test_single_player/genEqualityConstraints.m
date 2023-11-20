function constraint = genEqualityConstraints(Index)
    function h = constraintfunc(z)
%         h = []; %zeros(11*4,1);
        x = Index.getVar(join(['x',string(0)],''), 0, z);
        h = x;
        
        for t = 0:9
            x = Index.getVar(join(['x',string(t)],''), 0, z);
            u = Index.getVar(join(['u',string(t)],''), 0, z);
            xx = Index.getVar(join(['x',string(t+1)],''), 0, z);
            h(4*(t+1)+1:4*(t+1)+4) = xx - dynamics(x,u);
        end
    end
    constraint = @constraintfunc;
end


