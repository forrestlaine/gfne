function constraint = genInequalityConstraints(~)
    function g = constraintfunc(~)
        g = zeros(0,1);
    end
    constraint = @constraintfunc;
end