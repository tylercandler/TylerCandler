%% ASEN 4057 - HW 1 Q2- Main
%
% Author: Tyler Candler
% Collaborators: N/A
% Date: 1/12/2022

%% Housekeeping 

clear all; 
close all; 
clc;

%% Sample Data

%Constants
r = 3; %m
MW = 4.02;
W_payload = 5; %kg
W_dry = 0.6; %0.6

max_altitude = max_alt(r, W_payload, W_dry, MW);

fprintf("The max altitude the balloon can reach is %f meters", max_altitude)






%% Part 1 - Function

function weight = balloon_weight(radius, W_payload, W_dry_balloon, W_gas_molec)

%Inputs: Radius of balloon, weight of payload, dry weight of balloon, and
%the molecular weight of the gas in the balloon
%Outputs: Total weight of balloon
%Note: all weights are in kilograms, and all distances are in meters

rho_0 = 1.225; %kg/m^3 density at sea level
%Weight of the gas
W_gas = (4*pi*rho_0*radius^3*W_gas_molec)/(3*28.966);

%Summing weights
weight = W_gas + W_dry_balloon + W_payload;
end

%% Part 2 - Function


function weight = displacement_weight(radius, altitude)

%Inputs: Radius of balloon, current altitude of balloon
%Outputs: Total weight of the displaced gas
%Note: all weights are in kilograms, and all distances are in meters

%Solve for Temperature(C) and Pressure(kPa) based on current altitude
if 0 < altitude <= 11000
    T = 15.04 - 0.00649*altitude;
    P = 101.29 *((T + 273.1)/(288.08))^5.256;
elseif 11000 < altitude <= 25000
    T = -56.46;
    P = 22.65*e^(1.73-0.000157*h);
elseif altitude > 25000
    T = -131.21 + 0.00299*h;
    P = 2.488 * ((T + 273.1)/216.6)^-11.388;
end
%Solve for air density based on current Temp and Pressure
rho = P/(0.2869*(T + 273.1));


weight = (4*pi*rho*radius^3)/3;
end

%% Part 3 - Function

function h = max_alt(radius, W_payload, W_dry_baloon, W_gas_molec)
%Inputs: Radius of balloon, weight of payload, dry weight of balloon, and
%the molecular weight of the gas in the balloon
%Outputs: Max attainable height of balloon
%Note: all weights are in kilograms, and all distances are in meters

%Initial conditions at altitude of 0
h = 0;
W_air = displacement_weight(radius, h);
W_balloon = balloon_weight(radius, W_payload, W_dry_baloon, W_gas_molec);
    
    %Calculate whether balloon ways less than the displaced air
    while W_balloon < W_air
       h = h + 10;
       W_air = displacement_weight(radius, h);
    %  W_balloon = balloon_weight(radius, W_payload, W_dry_baloon, W_gas_molec);
    end

end



