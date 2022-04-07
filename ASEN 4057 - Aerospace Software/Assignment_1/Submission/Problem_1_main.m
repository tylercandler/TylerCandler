%% ASEN 4057 - HW 1 Q1 - Main
%
% Author: Tyler Candler
% Collaborators: N/A
% Date: 1/12/2022

%% Housekeeping 

clear all; 
close all; 
clc;

%% Constants

g = -9.81; % m/s

t = linspace (0,100, 1000); % time vector (arbitrary length)

V = input("What is the launch speed, in m/s? ");
theta = input("What is the launch angle, in degrees? ");
% V = 100;
% theta = 45;

%Calculate x position
x_vec = V * cosd(theta).*t;
%Calculate y position
y_vec = (0.5 * g * t.^2) + (V * sind(theta) * t);
%Find point at which y position becomes negative
y_neg = find(y_vec<0);


%Plot results
figure(1)

plot(x_vec(1: y_neg(1)-1),y_vec(1: y_neg(1)-1))
title("2D Position of a Projectile")
xlabel("X Position [m]")
ylabel("Y Position [m]")
yline(0)







