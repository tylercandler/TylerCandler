%% ASEN 4057 - HW 1 Q3 - Main
%
% Author: Tyler Candler
% Collaborators: N/A
% Date: 1/12/2022

%% Housekeeping 

clear all; 
close all; 
clc;

%% Sample Data
x = [-3, 12.1, 20, 0, 8, 3.7, -5.6, 0.5, 5.8, 10];
y = [-11.06, 58.95, 109.73, 3.15, 44.83, 21.29, -27.29, 5.11, 34.01, 43.25];
N = length(x);

%% Least Squares Estimate
A = sum(x);
B = sum(y);
C = sum(x.*y);
D = sum(x.*x);

m_sampledata = (A*B - N*C)/(A^2 - N*D)
b_sampledata = (A*C - B*D)/(A^2-N*D)

%% Lift Data
%clear least squares variables for next calculation
clear all;
%Experimental Data
alpha = [-5, -2, 0, 2, 3, 5, 7, 10, 14];
C_l = [-0.008, -0.003, 0.001, 0.005, 0.007, 0.006, 0.009, 0.017, 0.019];
N = length(alpha);

%Least Squares Calculation
A = sum(alpha);
B = sum(C_l);
C = sum(alpha.*C_l);
D = sum(alpha.*alpha);

m_lift = (A*B - N*C)/(A^2 - N*D)
b_lift = (A*C - B*D)/(A^2-N*D)

%Least Squares plottable line
x = linspace(-6,14,1000);
y = m_lift*x+b_lift;

%Plotting Results
figure(1)
hold on
scatter(alpha,C_l)
plot(x,y)
title("Coefficient of Lift vs Angle of Attack for F-117 Nighthawk ")
xlabel("Angle of Attack")
ylabel("Coefficient of Lift")
legend("Experimental Data","Least Squares Estimate",'Location','NW')


