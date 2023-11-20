% RLC Control
close all; clear; clc;

A = [0.983 0.1563; -0.1563 0.5921]; 
B = [0.0170; 0.1563];

rank([B A*B]);

T = 100;
x0 = [3; 1];

X = zeros(2,T+1);
X(:,1) = x0;

for t = 1:T
    X(:,t+1) = A*X(:,t);
end

time = linspace(0,(T+1)*0.2,T+1);

figure;
plot(time, X(1,:), 'linewidth', 3);
hold on;
plot(time, X(2,:), 'linewidth', 3);
xlabel('Time');
legend('Voltage over C','Current');

%%

X_ref = zeros(2,T+1);
X_ref(:,1) = x0;
for t = 1:T
    X_ref(1,t+1) = X_ref(1,1) + sin(t/5);
    X_ref(2,t+1) = X_ref(2,1) - sin(t/5);
end

figure;
hold on;
plot(time, X_ref(1,:),':b', 'linewidth', 3);
plot(time, X_ref(2,:),':r','linewidth', 3);
axis([0 T*0.2 -2 8])
xlabel('Time');
legend('Voltage ref','Current ref');

%% 

AA = zeros(T*2,T*3);
ind = 1;
for t = 1:T
    if (t>1)
        AA((t-1)*2+1:t*2,ind-2:ind-1) = A;
    end
    AA((t-1)*2+1:t*2,ind) = B;
    AA((t-1)*2+1:t*2,ind+1:ind+2) = -eye(2);
    ind = ind+3;
end

bb = zeros(T*2,1);
for t = 1:T
    bb(t*2-1:t*2) = X_ref(:,t+1)-A*X_ref(:,t);
end

%%
soln = AA'/(AA*AA')*bb;
X_soln = zeros(2,T+1);
U_soln = zeros(1,T);
X_soln(:,1) = x0;
ind = 1;
for t = 1:T
    U_soln(t) = soln(ind);
    ind = ind+1;
    X_soln(:,t+1) = soln(ind:ind+1)+X_ref(:,t+1);
    ind = ind+2;
end
%%
plot(time, X_soln(1,:),'b', 'linewidth', 2);
plot(time, X_soln(2,:),'r','linewidth', 2);

