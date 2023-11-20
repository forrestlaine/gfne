function xx = dynamics(x, u)
    xx = x + 0.1 * [x(3); x(4); u(1); u(2)];
end

