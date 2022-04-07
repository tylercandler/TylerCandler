%% ASEN 3111 - Computational Assignment 02
% % By Tyler Candler
% % Collaborators: Z. Vanlangendonck
% ASEN 3111 CA2 Main Script
% Feb 2022
clear all; close all; clc;
tic
%% Problem 1
disp("Starting Problem 1")
%% Flow Conditions
c = 2; %m chord length
alpha = 9; %degrees angle of attack
V_inf = 60; %m/s freestream velocity
p_inf = 85.5e3; %Pa freestream pressure
rho_inf = 1; %kg/m^3 air density

%% Compute coverged case for 1000 vortices 
    [x_converged, y_converged, u_converged, v_converged, p_converged, psi_converged, phi_converged, x_converged, y_converged] = PlotThinAirfoil(c, alpha, V_inf,  p_inf, rho_inf, 1000, 0);
    % Total velocity
    Vtot_converged = sqrt(u_converged.^2+v_converged.^2);
    % Pressure Coef
    Cp_converged = 1-(Vtot_converged/V_inf).^2;

%% For loop for 10 to 1000 vortices

N_vor = linspace(10,580,58); 
% I already calculated that the flow converges with 560 vortices,
% so I capped this vector at 580 to save time when running
for i = 1:length(N_vor)
    fprintf(['Calculating flow with ' num2str(N_vor(i)) ' out of ' num2str(N_vor(end)) ' vortices. \n'])
    [x, y, u, v, p, psi, phi, x_vortex, y_vortex] = PlotThinAirfoil(c, alpha, V_inf, p_inf, rho_inf, N_vor(i), 0);
   
    %Total Velocity
    Vtot = sqrt(u.^2+v.^2);
    % Pressure Coef from Vtot
    Cp = 1-(Vtot/V_inf).^2;
    
    %Compute root mean squared error values for V_inf, Cp, Psi, and Phi
    Vtot_RMS(i) = mean(mean(sqrt(((Vtot-Vtot_converged)./Vtot_converged).^2)));
    Cp_RMS(i) = mean(mean(sqrt(((Cp-Cp_converged)./Cp_converged).^2)));
    Psi_RMS(i) = mean(mean(sqrt(((psi-psi_converged)./psi_converged).^2)));
    Phi_RMS(i) = mean(mean(sqrt(((phi-phi_converged)./phi_converged).^2)));
end

%Calculate index where 1% error occurs for all stream variables
% Velocity 
Vtot_RMS_1 = lt(Vtot_RMS,0.01); % returns a logical array where the Vtot_RMS vec is less than 1%
i1_Vtot = find(Vtot_RMS_1, 1, 'first'); %first index that has a positive index, corresponds to the first index where error is less than 1%
N1_Vtot = N_vor(i1_Vtot); %use this index to find the number of vortices for this error
% Cp
Cp_RMS_1 = lt(Cp_RMS,0.01); %Repeat the steps above for all variables
i1_Cp = find(Cp_RMS_1, 1, 'first');
N1_Cp = N_vor(i1_Cp);
% Psi
Psi_RMS_1 = lt(Psi_RMS,0.01);
i1_Psi = find(Psi_RMS_1, 1, 'first');
N1_Psi = N_vor(i1_Psi);
% Phi
Phi_RMS_1 = lt(Phi_RMS,0.01);
i1_Phi = find(Phi_RMS_1, 1, 'first');
N1_Phi = N_vor(i1_Phi);


% Find the maximum number of vortexes for the above values that were
% calculated
N_max_vortex = max([N1_Phi, N1_Psi, N1_Cp, N1_Vtot]);
fprintf('The simulation converges to within 1 percent error using %i vortices. \n', N_max_vortex);

%Plot error
figure(1)
% figure('Units','inches','Position',[2 2 8 4]);
semilogy(N_vor, Cp_RMS, '-b');
hold on
semilogy(N_vor, Vtot_RMS, '-r');
semilogy(N_vor, Psi_RMS, '-g');
semilogy(N_vor, Phi_RMS, '-c');
semilogy([N_max_vortex N_max_vortex],[10^-6 10^0],'-k');
xlabel('Number of vortices');
ylabel('RMS Percent Error');
title('1 Percent Error Convergence');
grid on
% l = legend('C_p','U_{tot}','\Psi','\Phi');
% l.Interpreter = 'tex';
legend("Cp","Vtot","Psi","Phi",'Interpreter','tex')


%Plot streamlines, equipotential lines, and pressure contours
%run PlotThinAirfoil with flag set to 1 and N set to N_max_vortex
[x_final, y_final, u_final, v_final, p_final, psi_final, phi_final, x_vortices_final, y_vortices_final] = PlotThinAirfoil(c, alpha, V_inf,  p_inf, rho_inf, N_max_vortex, 1);
  
%% Problem 2
%Airflow over a thick airfoil using Vortex_Panel.m
clear all;
c = 2; %m chord length
alpha = 9; %degrees angle of attack
V_inf = 60; %m/s freestream velocity
p_inf = 85.5e3; %Pa freestream pressure
rho_inf = 1; %kg/m^3 air density

% Kuethe and Chowe Data

X_test = [1 .933 .750 .500 .250 .067 0 .067 .250 .500 .750 .933 1];
Y_test = [0 -.005 -.017 -.033 -.042 -.033 0 .045 .076 .072 .044 .013 0];

CP_KC = [0.2630 0.1969 0.2097 0.2667 0.4707 0.9929 -1.8101 -1.5088 -0.9334 -0.5099 -0.1688 0.1674];

%Verify Cp for Kuethe and Chowe Data

[Cl_test,Cp_test, scrap1, scrap2, scrap3] = Vortex_Panel(X_test, Y_test, V_inf, alpha, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% did not finish this problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Problem 3
disp("Starting Problem 3")
% Using the provided vortex panel method MATLAB function, preform a quantitative study
% of the convergence in the lift as a function of the panel resolution. 

clear all;
c = 2; %m chord length
alpha = 9; %degrees angle of attack
V_inf = 60; %m/s freestream velocity
p_inf = 85.5e3; %Pa freestream pressure
rho_inf = 1; %kg/m^3 air density


%% NACA 0012
m_0012 = 0/100; % maximum camber
p_0012 = 0/10; % location of maximum camber
t_0012 = 12/100; % thickness
error__0012 = 1; %initial error of 100 percent
i = 1;

N_0012 = linspace(10,1000);
% Cl and Cp with one Panel
[x_0,y_0] = NACA_SHAPE(m_0012, p_0012, t_0012, c, N_0012(i));
[Cl_0, Cp_0, scrap1, scrap2, scrap3] = Vortex_Panel(x_0, y_0, V_inf, alpha, 0); 
%while loop ends when error dips below 1 percent
%error is found by comparing the Cl value from the current iteration to
%that of the previous iteration
disp("Beginning Calculations")
while i ~= 101
    N = N_0012(i);
%     fprintf("Calculating %f panels \n",N_0012(i))
    [x_current,y_current] = NACA_SHAPE(m_0012, p_0012, t_0012, c, N); % (x,y) coordinates of the airfoil surface for current number of panels
    [Cl_current, Cp_current, scrap1, scrap2, scrap3] = Vortex_Panel(x_current, y_current, V_inf, alpha, 0); % Cl and Cp of the airfoil for current number of panels
    %store coefficient of lift for plotting purposes
    Cl_0012_vec(i) = Cl_current;

    %step once
    i = i+1;

end

%error less than 1 percent and associated index and number of panels
error_vec = abs(Cl_0012_vec-Cl_0012_vec(end))./Cl_0012_vec(end)*100;
error_vec = error_vec<0.01;
error_vec_index = find(error_vec,1,"first");
N_1percent = N_0012(error_vec_index);

%Display convergence to command prompt and plot Cl vs N
fprintf("Error converges to less than 1 percent using %f panels \n",N_1percent)
figure
plot(N_0012,Cl_0012_vec,'LineWidth',2)
hold on
xlabel("Number of Panels")
ylabel("Sectional Coefficient of Lift")
title("NACA-0012: Cl vs Number of Panels")
xline(N_1percent,'LineWidth',2)
legend("Data","Convergence Point",'Location','southeast')

%% Part 2
% Calculate Cl for a variety of alphas for 3 NACA airfoils and plot them using Vortex_Panel 
% Initialize Alpha Vector
alpha_vec = -5:5;
% N_1percent = 500; %temp
%% NACA 0006
m_0006 = 0/100; % maximum camber
p_0006 = 0/10; % location of maximum camber
t_0006 = 6/100; % thickness
[x_0006,y_0006] = NACA_SHAPE(m_0006, p_0006, t_0006, c, N_1percent);
for i = 1:length(alpha_vec)
    [Cl_0006(i), scrap4, scrap1, scrap2, scrap3] = Vortex_Panel(x_0006, y_0006, V_inf, alpha_vec(i), 0);
end

%% NACA 0024

[x_0012,y_0012] = NACA_SHAPE(m_0012, p_0012, t_0012, c, N_1percent);
for i = 1:length(alpha_vec)
    [Cl_0012(i), scrap4, scrap1, scrap2, scrap3] = Vortex_Panel(x_0012, y_0012, V_inf, alpha_vec(i), 0);
end


%% NACA 0024
m_0024 = 0/100; % maximum camber
p_0024 = 0/10; % location of maximum camber
t_0024 = 24/100; % thickness

[x_0024,y_0024] = NACA_SHAPE(m_0024, p_0024, t_0024, c, N_1percent);
for i = 1:length(alpha_vec)
    [Cl_0024(i), scrap4, scrap1, scrap2, scrap3] = Vortex_Panel(x_0024, y_0024, V_inf, alpha_vec(i), 0);
end

%plotting
figure
hold on
plot(alpha_vec,Cl_0006,'LineWidth',2)
plot(alpha_vec,Cl_0012,'LineWidth',2)
plot(alpha_vec,Cl_0024,'LineWidth',2)
title("Sectional Coeffiecient of Lift vs Angle of Attack for Several NACA Airfoils")
xlabel("Angle of Attack in degrees")
ylabel("Sectional Coefficient of Lift")
legend("NACA-0006","NACA-0012","NACA-0024")


%Calculate lift slope
p_0006 = polyfit(alpha_vec,Cl_0006,1)
p_0012 = polyfit(alpha_vec,Cl_0012,1)
p_0024 = polyfit(alpha_vec,Cl_0024,1)

fprintf("NACA-0006 Lift slop is approximately %f\n",p_0006)
fprintf("NACA-0012 Lift slop is approximately %f\n",p_0012)
fprintf("NACA-0024 Lift slop is approximately %f\n",p_0024)
fprintf("All airfoils have approximately zero lift at 0 degrees angle of attack \n")




toc
%% NACA__SHAPE Function Handle
function [xb, yb] = NACA_SHAPE(m,p,t,c,N)
% Inputs: m = max thickness
%         p = location of max thickenss
%         t = thickness
%         c = chord length
%         N = number of panels
% Outputs:
%         Boundary points along airfoil x and y vector
% Collaborators: Z. Vanlangendonck

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
end


%% PlotThickAirfoil Functiion Handle

function [] = PlotThickAirfoil(c, alpha, V_inf, p_inf, rho_inf, N, flag)



end
%% PlotThinAirfoil Function Handle

function [x, y, u, v, p, psi, phi, x_vortices, y_vortices] = PlotThinAirfoil(c, alpha, V_inf, p_inf, rho_inf, N, flag)
% Inputs: flow condiitions: c, alpha, V_inf, p_inf, rho_inf
%                           N: number of vortices
%                           flag: 1 to plot, 0 to not plot
% Outputs: stream function, locations of vortices, plots

%Collaborators: N/A, but heavy influence from the Lifting_Cylinder.m script
%that was provided



resolution = 0.02;
%resolution for meshgrid

%% Divide chord into separate vortices
dx = c/N; %distance between vortices
x_vortices = linspace(0.5*dx,c-0.5*dx,N); % spacing vortices along chord, starting and ending at 0.5dx from the LE and TE, with N points
y_vortices = zeros(size(x_vortices));
Gamma = 2*alpha*V_inf*sqrt((1-x_vortices/c)./(x_vortices/c))*dx;

%% Define Domain
xmin=-1*c;
xmax=3*c;
ymin=-1*c;
ymax=1*c;

%% Define Number of Grid Points
nx=(xmax-xmin)/resolution; % steps in the x direction
ny=(ymax-ymin)/resolution; % steps in the y direction

%% Create mesh over domain using number of grid points specified
[x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));


%% Add the Uniform flow contribution
u = V_inf*cos(alpha);
v = V_inf*sin(alpha);
psi = u.*y-v.*x;
phi = u.*x+v.*y;

%% for loop over all the vortices
for  i = 1:N
    %compute influence of each individual vortex
    
    %radius and angle to each point
    radius = sqrt((x-x_vortices(i)).^2+(y-y_vortices(i)).^2);
%   theta = floor(atan2((y-y_vortices(i)),(x-x_vortices(i)))+alpha,pi*2);
    theta = mod(atan2((y-y_vortices(i)),(x-x_vortices(i)))+alpha,pi*2);
    %Stream function, velocity components
    phi_vortex = - Gamma(i) * theta/(2*pi);
    psi_vortex = Gamma(i)*log(radius)/(2*pi);
    u_t = -Gamma(i)./(2*pi*radius);
    u_vortex = -u_t.*sin(theta);
    v_vortex = u_t.*cos(theta);

    % Sum all of these values
    phi = phi+phi_vortex;
    psi = psi+psi_vortex;
    u = u+u_vortex;
    v = v+v_vortex;

end
   p = p_inf;   
%if flag delimmeter is 1, plot equipotential lines, pressure contours and streamlines 
if flag == 1
    
    %Total Velocity
    Vtot = sqrt(u.^2+v.^2);
    % Pressure Coef from Vtot
    Cp = 1-(Vtot/V_inf).^2;
        
    %Plot Streamlines
    figure
    hold on
    plot(x_vortices,y_vortices,'LineWidth',2);
    contour(x,y,psi,50)
    axis([-2 4 -2 2]);
    axis equal
    xlabel('x, [m]');
    ylabel('y, [m]');
    title(['Stream Function with ' num2str(N) ' Vortices'])


    figure
    % Plot equipotential lines
    
    hold on
    plot(x_vortices,y_vortices,'LineWidth',2);
    contour(x,y,phi, 50)
    axis([-2 4 -2 2]);
    axis equal
    xlabel('x, [m]');
    ylabel('y, [m]');
    title(['Equipotential Lines with ' num2str(N) ' Vortices']);

    % Plot Cp
    figure
    hold on
    plot(x_vortices,y_vortices,'LineWidth',2);
    contour(x,y,Cp,(linspace(-2,1,50)));
    axis([-2 4 -2 2]);
    axis equal
    xlabel('x, [m]');
    ylabel('y, [m]');
    title(['Pressure Coefficient Contours with ' num2str(N) ' Vortices']);
    cb = colorbar;
    cb.Label.String = 'C_p';


    end
end