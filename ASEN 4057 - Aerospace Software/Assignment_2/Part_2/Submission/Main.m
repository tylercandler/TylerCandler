%% ASEN 4057 Homework #2 
%Tyler Candler & Hannah Johnson

profile on
n = 100;
M = magic(n);

p = profile('info')
save myprofiledata p
%% Housekeeping
clear all; close all; clc;

%% Constants

G=6.67*10^(-11);
mM = 7.34767309*10^22; % mass of moon in kg
mE = 5.97219*10^24; % mass of Earth kg
mS = 28833; % mass of spacecraft in kg
rM = 1737100; %radius of moon in m
rE = 6371000; %radius of the Earth in m
dEM=384403000; %distance from Earth to moon in m




%% Initial Conditions
dES_0 = 340000000; 
vS_0 = 1000;
theta_s = 50;
xS_0 = dES_0 * cosd(theta_s); 
yS_0 = dES_0 * sind(theta_s); 
vSx_0 = vS_0 * cosd(theta_s); 
vSy_0 = vS_0 * sind(theta_s); 



dEM_0 = 384403000; 
vM_0 = sqrt((G*mE^2)/((mE+mM)*dEM_0));
theta_m = 42.5;
xM_0 = dEM_0 * cosd(theta_m); 
yM_0 = dEM_0 * sind(theta_m); 
vMx_0 = -vM_0 * sind(theta_m); 
vMy_0 = vM_0 * cosd(theta_m); 

xE_0 = 0;  
yE_0 = 0;  
vEx_0 = 0;  
vEy_0 = 0; 


% IC = [xS_0, yS_0, vSx_0, vSy_0, xM_0, yM_0, vMx_0, vMy_0, xE_0, yE_0, vEx_0, vEy_0]';%vector of initial conditions

options = odeset('Events',@myevents,'Reltol',1e-15,'AbsTol',1e-15); %options to set terminal conditions


%% Optimization
opt2 = optimset('TolFun', 1e-4, 'TolX', 1e-5, 'Display','iter','MaxFunEvals',1e5); %options for minimizing function

% Optimize to find minimum delta V
delta_V_S_0 = [1000,1000]; %Initial guess, arbitrary but the closer it is, the faster the grid search is function works


%Grid search to further narrow down the initial guess (step size of 10)
max_delta_V_S = 100; % max delta V as per constraints
%Upper hemisphere
for delta_V_X=-100:20:100
    for delta_V_Y =0:20:100
        if sqrt(delta_V_X^2 + delta_V_Y^2)<max_delta_V_S
            v = OptimizationFunction([delta_V_X, delta_V_Y]);
            if v < 100
                max_delta_V_S = v;
            end
            if v<norm(delta_V_S_0)
                delta_V_S_0 = [delta_V_X, delta_V_Y];
            end
        end
    end
end
% Lower hemisphere
for delta_V_X=-100:20:100
    for delta_V_Y= 0:-1*20:-100
        if sqrt(delta_V_X^2 + delta_V_Y^2)<max_delta_V_S
            v = OptimizationFunction([delta_V_X, delta_V_Y]);
            if v < 100
                max_delta_V_S = v;
            end
            if v<norm(delta_V_S_0)
                delta_V_S_0=[delta_V_X, delta_V_Y];
            end
        end
    end
end
% delta_V_S_0

%Use fminsearch to precisely find the minimum delta_V
delta_V_S = fminsearch(@(delta_V_S)OptimizationFunction(delta_V_S),delta_V_S_0, opt2);







%%  ODE call with added delta V
%tspan, arbitrary time as the ode should end based on end conditions
tspan = [0 10000000];

% IC = [xS_0, yS_0, vSx_0 + delta_V_S(1), vSy_0 + 49, xM_0, yM_0, vMx_0, vMy_0, xE_0, yE_0, vEx_0, vEy_0]';%vector of initial conditions

IC_new = [xS_0, yS_0, vSx_0 + delta_V_S(1), vSy_0 + delta_V_S(2), xM_0, yM_0, vMx_0, vMy_0, xE_0, yE_0, vEx_0, vEy_0]';%vector of initial conditions


ODEoptions = odeset('Events', @myevents, 'RelTol', 1e-8, 'MaxStep', 1e5);
% ode call to integrate our Equations of Motion
[t,out] = ode45(@(t,input) ODEFUN(t,input),tspan,IC_new,ODEoptions);




%% Output vector

%Extract positions of the moon, satellite, and Earth from the ODE call

xS = (out(:,1));
yS = (out(:,2));
xM = out(:,5);
yM = out(:,6);
xE = out(:,9);
yE = out(:,10);


%% Plotting
figure(1)
hold on %I'm trying dear god I am trying

%plot spaceship's path
plot(xS,yS,'r')

%plot path of the moon
plot(xM,yM,'b')



%Plot Earth
%used given radius and assumed that the Earth is relatively stationary
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', rE);
plot(p, 'FaceColor', [0.5 0.7 0.8])
axis equal

title('Trajectory more like tragedy')
xlabel('Meters');
ylabel('Meters')
legend('Spaceship Trajectory','Moon Trajectory','Earth','location','southeast')
fprintf("The change in velocity required to get the spaceship back to Earth is %f m/s in the x direction and %f m/s in the y direction",delta_V_S(1),delta_V_S(2))


%% Optimization Function
function result = OptimizationFunction(delta_V_S)
% Constants
G=6.67*10^(-11);
mM = 7.34767309*10^22; % mass of moon in kg
mE = 5.97219*10^24; % mass of Earth kg
mS = 28833; % mass of spacecraft in kg
rM = 1737100; %radius of moon in m
rE = 6371000; %radius of the Earth in m
dEM=384403000; %distance from Earth to moon in m

% Initial Conditions

%Reset initial conditions
dES_0 = 340000000; 
vS_0 = 1000;
theta_s = 50;
xS_0 = dES_0 * cosd(theta_s); 
yS_0 = dES_0 * sind(theta_s); 
vSx_0 = vS_0 * cosd(theta_s); 
vSy_0 = vS_0 * sind(theta_s); 



dEM_0 = 384403000; 
vM_0 = sqrt((G*mE^2)/((mE+mM)*dEM_0));
theta_m = 42.5;
xM_0 = dEM_0 * cosd(theta_m); 
yM_0 = dEM_0 * sind(theta_m); 
vMx_0 = -vM_0 * sind(theta_m); 
vMy_0 = vM_0 * cosd(theta_m); 

xE_0 = 0;  
yE_0 = 0;  
vEx_0 = 0;  
vEy_0 = 0; 

% Declare initial condition vector with added delta V
IC = [xS_0, yS_0, vSx_0 + delta_V_S(1), vSy_0+ delta_V_S(2), xM_0, yM_0, vMx_0, vMy_0, xE_0, yE_0, vEx_0, vEy_0]';%vector of initial conditions

% ODE CALL
options = odeset('Events', @myevents, 'RelTol', 1e-8, 'MaxStep', 1e5);  % custom options
tspan = [0, 10^10];       %integration tspan
[t, out, te,ye,ie] = ode45(@(t,input)ODEFUN(t, input),tspan, IC, options); 
% [t,out] = ode45(@(t,input) ODEFUN(t,input),tspan,IC,options);

% if simulation results in the spaceship returning back to earth, return
% the total delta V that resulted in that successful simulation, otherwise
% return very high value
if ie==2
    result = sqrt(delta_V_S(1)^2 + delta_V_S(2)^2);  % Magnitude of delta_V_S vector
else
    result=100000;
end


if ie==1
    fprintf('Spacecraft hit the moon \n');
end
if ie==2
    fprintf('Spacecraft returned to Earth \n');
end
if ie==3
    fprintf('Spaceship was lost \n');
end
end




%% End conditions function
function [value,isterminal,direction]=myevents(t,y)
rE=6371000;
rM=1737100;

%distance between spaceship and moon
dSM=sqrt((y(1)-y(5))^2+(y(2)-y(6))^2);

%distance between Earth and spaceship
dES=sqrt((y(1)-y(9))^2+(y(2)-y(10))^2); 

%distance between Earth and Moon
dEM=sqrt((y(9)-y(5))^2+(y(10)-y(6))^2);

% Event 1: Spacecraft crashes into moon
value(1)      = 1*(dSM>rM);
isterminal(1) = 1;
direction(1)  = 0;


% Event 2: Spacecraft returns to Earth
value(2)      = 1*(dES>rE);    
isterminal(2) = 1;
direction(2)  = 0;


% Event 3: Spacecraft is Lost
value(3)      = 1*(dES < 2*dEM);
isterminal(3) = 1;
direction(3)  = 0;

end
