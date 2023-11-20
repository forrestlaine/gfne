clear; clc; close all;

T = readtable('~/Downloads/Midterm_1_PDF_scores.csv');

scores = T.TotalScore;
scores = scores(~isnan(scores));

figure;
hist(scores,30);

figure;
hist(scores+22,30);