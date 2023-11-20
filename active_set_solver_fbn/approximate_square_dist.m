% test approximate square dist

w1 = 3;
h1 = 1;

w2 = 5;
h2 = 3;

x1 = [0;0];
x2 = [4;3];

delta = x1-x2;
dist = sqrt(delta'*delta);

d1 = 0;
d2 = 0;
if abs(delta(1)) > w1
    d1 = w1 / abs(delta(1)) * dist;
elseif abs(delta(2)) > h1
    d1 = h1 / abs(delta(2)) * dist;
else
    d1 = dist;
end
if abs(delta(1)) > w2
    d2 = w2 / abs(delta(1)) * dist;
elseif abs(delta(2)) > h2
    d2 = h2 / abs(delta(2)) * dist;
else
    d2 = dist;
end

approximate_dist = max(0,dist-d1-d2);

H = [eye(2) -eye(2);
    -eye(2) eye(2)];
f = [0;0;0;0];

lb = [x1(1)-w1; x1(2)-h1; x2(1)-w2; x2(2)-h2];
ub = [x1(1)+w1; x1(2)+h1; x2(1)+w2; x2(2)+h2];

[z,fval] = quadprog(H,f,[],[],[],[],lb,ub);
true_dist = sqrt(2*fval);



