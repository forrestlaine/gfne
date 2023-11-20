function [c, ceq, gradc, gradceq] = con_helper_eval(x, ch, ceqh, gradch, gradceqh)
    c = sparse(ch(x));
    ceq = sparse(ceqh(x));
    if nargout > 2
        gradc = sparse(gradch(x)');
        gradceq = sparse(gradceqh(x)');
    end
end