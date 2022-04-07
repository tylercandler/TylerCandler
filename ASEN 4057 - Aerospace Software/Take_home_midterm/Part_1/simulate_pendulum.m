%% Function Header
function [t,y] = simulate_pendulum(tspan, y0, l1, l2, m1, m2)
% y0 must be of the form IC = [theta1_0, v_theta1_0, theta2_0, v_theta2_0]'
L1 = l1; % lazy
L2 = l2; % lazy
g = 9.81;
%% ODE45 Call
[t, y] = ode45(@(t,input) EOM(t,g,L1,L2,m1,m2,input),tspan,y0);
end