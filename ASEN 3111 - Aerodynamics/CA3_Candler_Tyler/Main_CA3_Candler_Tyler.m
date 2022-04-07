%% ASEN 3111 - Computational Assignment 3 - Main
% Employs thin airfoil theory to find coeffients of lift, as well as zero
% lift angle of attack and lift slope for different NACA airfoils
% Experiments with Prandtl's Lifting Line Theory to explore the effect of
% Aspect Ratio and Taper Ratio on Drag coefficient and span efficiency
% factor
%
% Author: Tyler Candler
% Collaborators: Z. Vanlangendonck, A. Gillepsie
% Date: 04/03/2022

close all; clear all; clc;

%% Problem 1
disp("Starting Problem 1")
%Declare NACA numbers as strings for later computation
NACA_0012 = '0012';
NACA_2412 = '2412';
NACA_4412 = '4412';

%Use NACA_Airfoil to find airfoil boundary conditions, zero lift AOA, and
%lift slope via thin airfoil theory
[xb1, yb1,aL0_1, m_cL] = NACA_Airfoil(NACA_0012, 2, 67);
[xb2, yb2,aL0_2, m_cL] = NACA_Airfoil(NACA_2412, 2, 67);
[xb3, yb3,aL0_3, m_cL] = NACA_Airfoil(NACA_4412, 2, 67);

%Use Vortex Panel Method to solve for zero lift AOA and lift slope for the
%same airfoils
[a_0_0012, a_L0_0012, b_0012] = LiftSlope(NACA_0012,2,100);
[a_0_2412, a_L0_2412, b_2412] = LiftSlope(NACA_2412,2,100);
[a_0_4412, a_L0_4412, b_4412] = LiftSlope(NACA_4412,2,100);


%Display results in a table
RESULT = ["Lift Slope (Thin Airfoil Theory)[1/rad]";"Zero Lift AOA (Thin Airfoil Theory)[deg]";"Lift Slope (Vortex Panel)[1/rad]";"Zero Lift AOA (Vortex Panel)[deg]"];
NACA_0012 = [m_cL, aL0_1*180/pi ,a_0_0012, a_L0_0012]';
NACA_2412 = [m_cL, aL0_2*180/pi ,a_0_2412, a_L0_2412]';
NACA_4412 = [m_cL, aL0_3*180/pi ,a_0_4412, a_L0_4412]';

data = table(RESULT, NACA_0012, NACA_2412, NACA_4412);
disp(data)

%Plot results
figure(1)
hold on
AOA = deg2rad(-5:10);


y_0012 = a_0_0012*AOA+b_0012;
plot(AOA, y_0012)

y_2412 = a_0_2412*AOA+b_2412;
plot(AOA, y_2412)

y_4412 = a_0_4412*AOA+b_4412;
plot(AOA, y_4412)
yline(0,"--")
xline(deg2rad(a_L0_0012),"--")
xline(deg2rad(a_L0_2412),"--")
xline(deg2rad(a_L0_4412),"--")
title("Coeffient of Lift vs AOA")
xlabel("AOA [rad]")
ylabel("Sectional Coefficient of Lift")
legend("NACA 0012","NACA 2412","NACA 4412","Zero Lift Angle of attack",'Location',"NW")



%% Problem 2
disp("Starting Problem 2")

%Givens for aircraft wing
b = 33+4/12; %span [ft]
c_r = 5+4/12; %root chord [ft]
c_t = 3+8.5/12; %tip chord [ft]
NACA_r = '2412'; %root airfoil
NACA_t = '0012'; %tip airfoil
geo_r = deg2rad(1); %root geometric angle of attack [rad]
geo_t = 0; %root geometric angle of attack [rad]
V_inf = 82*1.68781; %free-stream velocity [ft/s]
% rho_inf = 0.0023769;%free-stream air density [slugs/ft^3] at 10,000 ft alt
rho_inf = 17.56; %slugs
S = (1/2)*(c_t+c_r)*b; %surface area [ft^2]

% Calculate lift slope and zero-lift angle of attack of each NACA airfoil
[a0_r, aero_r] = LiftSlope(NACA_r,c_r,100); %root airfoil
[a0_t, aero_t] = LiftSlope(NACA_t,c_t,100); %tip airfoil

%convert to degrees
aero_r = deg2rad(aero_r);
aero_t = deg2rad(aero_t);

%Run PLLT with large number of panels to obtain converged solution
[e,C_L,C_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,5000);

L = (1/2)*C_L*rho_inf*S*V_inf^2;
D = (1/2)*C_Di*rho_inf*S*V_inf^2;
fprintf("Converged solution for the given airfoil: \ne = %f \nC_L = %f \nLift = %f pounds \nC_Di = %f \nDrag = %f pounds \n", e,C_L,L,C_Di,D)

%Error calculations

% Panels required for 10,1, and 0.1 percent error
% initialize vectors
N_vec = 2:10000; 
err = 100; 
c_L_vec = NaN(2000,1);
c_D_vec = NaN(2000,1);
err_vec = NaN(2000,1);
%for loop to find errors
for i = 1:length(N_vec)
    [scrap, c_L_vec(i), c_D_vec(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N_vec(i));
    err_vec(i) = 100*abs(C_L - c_L_vec(i))/C_L;
    if err_vec(i) < 0.005
        break
    end
end
    
c_L_vec = c_L_vec(~isnan(c_L_vec));
c_D_vec = c_D_vec(~isnan(c_D_vec));
err_vec = err_vec(~isnan(err_vec));
N_vec = N_vec(1:length(err_vec));
    
%10 percent error
index_10 = find(err_vec <= 10, 1);
%1 percent error
index_1 = find(err_vec <= 1, 1);
%Point 1 percent error
index_pt1 = find(err_vec <= 0.1, 1);
%Print results to command window
fprintf("Required number of odd terms for: \nTen percent error: %d panels \nOne percent error: %d panels \nPoint one percent error: %d panels \n", N_vec(index_10),N_vec(index_1),N_vec(index_pt1))


%% Problem 3
disp("Starting Problem 3")

AR_vec = [4,6,8,10];
taper_ratio = linspace(0,1,101);
e_mat = NaN(length(taper_ratio),length(AR_vec));
c_t_temp = c_t;
b_temp = b;

%         c_t = c_r*j;
%         b = i*((1+j)*c_r)/2;
%         [e_vals(k,l),~,~] = PLLT(b,2*pi,2*pi,c_t,c_r,0,0,1,1,N);
%         l = l + 1;

for Aspect = 1:length(AR_vec) %AR_vec
    for TR = 1:length(taper_ratio) %taper_ratio
        c_t_temp = c_r * taper_ratio(TR);
        b_temp = AR_vec(Aspect)*((1+taper_ratio(TR))*c_r)/2;
        [e_mat(TR,Aspect),~,~] = PLLT(b_temp,2*pi,2*pi,c_t_temp,c_r,0,0,1,1,20);       
    end
end

figure(2)
hold on;
plot(taper_ratio,e_mat(:,1))
plot(taper_ratio,e_mat(:,2))
plot(taper_ratio,e_mat(:,3))
plot(taper_ratio,e_mat(:,4))
title("Span Efficiency Factor vs Taper Ratio")
xlabel("Taper ratio, c_t/c_r")
ylabel("Span efficiency factor, e")
legend("AR = 4","AR = 6","AR = 8","AR = 10", 'Location', "SE")





%% NACA_Airfoil Function Handle
function [xb, yb, aL0, dCLda] = NACA_Airfoil(NACA_NUM,c,N)
% Inputs: m = max thickness
%         NACA_NUM =  4 character NACA number as a string
%                from which we can obtain:
%                m = max camber
%                p = location of max camber
%                t = thickness
%         c = chord length 
%         N = number of panels
% Outputs:
%         Boundary points along airfoil x and y vector
%         aL0 = zero lift angle of attack
%         dCLda = lift slope according to thin airfoil theory (2pi)
% Author: Tyler Candler
% Collaborators: Z. Vanlangendonck, A. Gillepsie
% Date: 04/04/2022

%Extract data from NACA number
m = str2num(NACA_NUM(1))/100;
p = str2num(NACA_NUM(2))/10;
t = str2num(NACA_NUM(3:4))/100;

%chord vector
x_chord = linspace(0,c,(N/2)+1);

% EQ from problem 2
y_calc = t/0.2*(0.2969*sqrt(x_chord/c)-0.1260*(x_chord/c)-0.3516*(x_chord/c).^2+0.2843*(x_chord/c).^3-0.1036*(x_chord/c).^4);

% Mean Camber Line
if m ~= 0 && p ~= 0 % if m and p are zero
    y_camber = zeros(1,length(x_chord)); %initialize y vector
    for i = 1:length(x_chord) % loop over whole x vector
        if x_chord(i) <= p*c % if x is before max camber
            y_camber(i) = m*(x_chord(i)/(p^2))*(2*p-(x_chord(i)/c));
        elseif x_chord(i) >= p*c % if x is after max camber
            y_camber(i) = m*((c-x_chord(i))/((1-p)^2))*(1+(x_chord(i)/c)-2*p);
        end
    end
else
    y_camber = zeros(1,length(x_chord));
end

% Zeta
zeta = atan(diff(y_camber));
zeta = [zeta, 0];

% Upper Surface
x_upper_bounds = x_chord - (y_calc.*sin(zeta));
y_uupper_bounds = y_camber + (y_calc.*cos(zeta));

% Lower Surface
x_lower_bounds = x_chord + (y_calc.*sin(zeta));
y_lower_bounds = y_camber - (y_calc.*cos(zeta));

% Add lower and upper surfaces in a counter clockwise direction
xb = [flip(x_lower_bounds), x_upper_bounds(2:end)];
yb = [flip(y_lower_bounds), y_uupper_bounds(2:end)];

%Lift slope
dCLda = 2*pi; %Thin Airfoil Theory


if m ~= 0 && p ~= 0 % if m and p are zero

    %Zero lift angle of attack
    % 

    theta_p = cos(1-2*p/c)^-1;

    % x = (c/2*(1-cos(theta))
    %Calculate dz/dx
    % dzdz1 = -m*x*(2*p-(x/c))/p^2;%before max camber
    %sub in theta
    % dzdz1 = -m*(c/2*(1-cos(theta)))*(2*p-(c/2*(1-cos(theta))/c))/p^2;
    fun1 = @(theta) -m.*(c/2.*(1-cos(theta))).*(2.*p-(c./2.*(1-cos(theta))/c))/p.^2.*(cos(theta)-1);

    % dzdx2 = -(2*m*(x-c*p))/(c*(p-1)^2);%after max camber
    %sub in theta
    % dzdx2 = -(2*m*((c/2*(1-cos(theta))-c*p)))/(c*(p-1)^2);%after max camber
    fun2 = @(theta) -(2.*m.*((c./2.*(1-cos(theta))-c.*p)))./(c.*(p-1).^2).*(cos(theta)-1);

    %integrate to calculate zero lift angle of attack

    integral1 = integral(fun1,0,theta_p);
    integral2 = integral(fun2,theta_p,pi);

    aL0 = -1/pi*(integral1 +integral2);
else
    aL0 = 0;
end

end
%% Lift Slope and Zero lift AOA function
function [a, aL0, b] = LiftSlope(NACA, c ,N)
%Calculates zero lift AOA and lift slope of a defined NACA airfoil using
%the vortex panel method
%Inputs:    NACA = NACA number as a string
%           c = chord length
%           N = number of panels to model airfoil
%Outputs:   a = lift slope
%           aL0 = zero lift AOA'
%           b = y intercept of the linear regression
% Author: Tyler Candler
% Collaborators: None
% Date: 04/04/2022

% Call NACA_Airfoil to find boundary points
[xb, yb, ~] = NACA_Airfoil(NACA,c,N);

% Cl at varying AOA
AOA_vec = deg2rad(-5:10); % AOA vec
Cl_vec = zeros(length(AOA_vec),1); % initialize Cl vec

for i = 1:length(AOA_vec)
    [Cl_vec(i),~] = Vortex_Panel(xb,yb,60,AOA_vec(i),0); % sectional Cl using Vortex_Panel
end
    
% Linear regression of lift slope
a_val = polyfit(AOA_vec,Cl_vec,1); % Linear regression using polyfit
a = a_val(1) * 60; % lift slope [1/rad]
aL0 = -a_val(2)/a_val(1); % zero-lift angle of attack [rad]
b = a_val(2);




end
%% Prandtl Lifting Line Theory Function

function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
% This function solves the Prandtly Lifitng Line Equation for finite wings
% w/ thick airfoils, using the vortex panel method to approximate lift
% Inputs:
%           b       - span
%           a0_t    - cross-sectional lift slope at the tip [1/rad]
%           a0_r    - cross-sectional lift slope at the root [1/rad]
%           c_t     - chord length at the tip
%           c_r     - chord length at the root
%           aero_t  - zero-lift angle of attack at the tip [rad]
%           aero_r  - zero-lift angle of attack at the root [rad]
%           geo_t   - geometric angle of attack at the tip [rad]
%           geo_r   - geometric angle of attack at the root [rad]
%           N       - number of odd terms to include in the series 
%                       expansion for circulation
% Outputs:
%           e       - span efficiency factor
%           c_L     - coefficient of lift
%           c_Di    - coefficient of induced drag
%
% Author: Tyler Candler
% Collaborators: None
% Date: 04/04/2022


%aspect ratio
S = (1/2)*(c_t+c_r)*b; %Planform area of the wing [ft^2]
AR = S^-1* b^2; %Aspect Ratio

%Theta and y
theta = (1:N)*pi/(2*N); %i*pi/2N
y = -b/2*cos(theta);

% Spanwise chord length, lift slope, and AOA
a0 = a0_r + (a0_t-a0_r)*y/(-b/2);
c = c_r + (c_t-c_r)*y/(-b/2);
aero = aero_r + (aero_t-aero_r)*y/(-b/2);
geo = geo_r + (geo_t-geo_r)*y/(-b/2);

%Solve for LHS and RHS
%Initialize arrays
LHS = NaN(N,N);
RHS = NaN(N,1);

for i = 1:N
    for j = 1:N
        LHS(i,j) = ((4*b)/(a0(i)*c(i)))*sin((2*j-1)*theta(i)) + (2*j-1)*sin((2*j-1)*theta(i))/sin(theta(i));
    end
    RHS(i) = geo(i) - aero(i);
end

%Solve for A coefficients
Aodd = LHS\RHS;

%Solve for Cl using coefficients
c_L = Aodd(1)*pi*AR;

%Solve for span efficiency factor
delta = 0;
for j = 2:N
    delta = delta + (2*j-1)*(Aodd(j)/Aodd(1))^2;
end
e = 1/(1+delta);

%Finally solve for induced drag coefficient
c_Di = (c_L)^2/(pi*e*AR);


end

