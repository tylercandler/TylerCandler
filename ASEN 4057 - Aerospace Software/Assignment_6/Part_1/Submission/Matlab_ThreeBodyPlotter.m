% Matlab_ThreeBodyPlotter.m
% ASEN 4057 Assigment 6
% Author: Dante Vigil

% Purpose: A simple plotting script of the ThreeBody_data.csv output from
% ThreeBody.c
clc, clear all, close all;
%% Data Load In
data = readmatrix("ThreeBody_data.csv");

dSM = sqrt((data(:,1)-data(:,5)).^2 + (data(:,2)-data(:,6)).^2);

figure
grid on
hold on

title("Trajectory")
xlabel("X [m]")
ylabel("Y [m]")
plot(data(:,1),data(:,2),'r--')
plot(data(:,5),data(:,6),'b--')

xlim([2*10^8 3.2*10^8])
ylim([2.3*10^8 3.2*10^8])

hold off