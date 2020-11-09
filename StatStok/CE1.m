%% Stationary stochastic processes, computer exercise 1
close all
clear all
clc
addpath('all_files')

% 2.1
disp('2.1')
% Load data
load data1.mat
% Plot
figure(2)
plot(data1.x)
hold on
plot(zeros(size(data1.x)),'r')
hold off
xlabel('Measurements')
title('2.1')
% Estimated expected value
m = mean(data1.x)
% Confidence interval
lower = m - 1.96*std(data1.x)/sqrt(length(data1.x));
upper = m + 1.96*std(data1.x)/sqrt(length(data1.x));
conf_int = [lower upper]

% 2.2
% Load and decide max time shift
load covProc
n = 3;
% Plot y(t) against y(t-k)
for k = 1:n
    figure(k+2)
    plot(covProc(1:end-k), covProc(k+1:end),'.')
    title(['2.2 Plot for k = ' num2str(k)])
    xlabel('y(t)')
    ylabel(['y(t-' num2str(k) ')'])
end
% Look at covariance function
[ycov,lags] = xcov(covProc,20,'biased');
figure(6)
subplot(2,1,1)
plot(lags,ycov,'o-')
hold on
plot(lags,zeros(size(ycov)),'r')
hold off
xlabel('Lags')
title('2.2 Covariance')
% Look at correlation function
[ycor,lags] = xcov(covProc,20,'coeff');
subplot(2,1,2)
plot(lags,ycor,'o-')
hold on
plot(lags,zeros(size(ycor)),'r')
hold off
xlabel('Lags')
title('2.2 Correlation')

% 2.3
disp('2.3')
% Create data
figure(1)
simuleraenkelsumma
% Open interface
%spekgui

% 3.1
% Load data
load cello.mat
load trombone.mat
% Open interface
%spekgui

% 3.2
% Down-sample cello
n = 2;
cello2.x = cello.x(1:n:end);
cello2.dt = cello.dt*n;
% Down-sample trombone
k = 4;
tromboneD.x = trombone.x(1:k:end);
tromboneD.dt = trombone.dt*k;
% Low-pass filter then down-sample
cello2_filt.x = decimate(cello.x,2);
cello2_filt.dt = cello.dt*2;
% Open interface
%spekgui


