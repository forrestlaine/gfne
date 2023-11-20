% Test measurement model
N = 1000;

s1 = 1;
m1 = -1;

s2 = .25;
s3 = .25;

u2 = 0.5;
u3 = 0.1;

a1 = -2.5;
a2 = -1.5;
a3 = -0.5;


X = linspace(-5,1,N);

for n = 1:N
    g1(n) = 1/sqrt(2*pi*s1*s1)*exp(-0.5*((X(n)-m1)/s1)^2);
    if X(n) < a1 
        g2(n) = 1 - (0.9*(1-u2)*exp(-0.5*((X(n)-a1)/s2)^2) + u2/2);
    elseif X(n) > a2
        g2(n) = 1 - (0.9*(1-u2)*exp(-0.5*((X(n)-a2)/s2)^2) + u2/2);
    else
        g2(n) = 1 - (0.9*(1-u2) + u2/2);
    end
    
    if X(n) < a2 
        g3(n) = 0.9*(1-u3)*exp(-0.5*((X(n)-a2)/s3)^2) + u3/2;
    elseif X(n) > a3
        g3(n) = 0.9*(1-u3)*exp(-0.5*((X(n)-a3)/s3)^2) + u3/2;
    else
        g3(n) = 0.9*(1-u3) + u3/2;
    end
    
    g4(n) = g1(n)*g2(n)*g3(n);
end

figure; hold on;
plot(X,g1);
plot(X,g2);
plot(X,g3);
plot(X,g4);