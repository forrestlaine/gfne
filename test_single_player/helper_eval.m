function [f, grad] = helper_eval(x, fh, gradh)
    f = full(fh(x));
    if nargout > 1
        grad = sparse(gradh(x));
    end
end