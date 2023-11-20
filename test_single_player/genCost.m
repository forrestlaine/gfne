function cost = genCost(Index)
    function out = costfunc(z)
        out = 0.0;
        for t = 0:9
            x = Index.getVar(join(['x',string(t)],''),0, z);
            u = Index.getVar(join(['u',string(t)],''),0, z);
            out = out + norm(x-[3,4,0,0]) + norm(u);
        end
        x = Index.getVar(join(['x',string(10)],''), 0, z);
        out = out + norm(x-[3,4,0,0]);
    end
    cost = @costfunc;
end

