clc; clear; close all;
X = linspace(-5,5,100);

y1 = X.^2 + 1;

figure;
hold on;
plot(X,y1);
plot(X,X);