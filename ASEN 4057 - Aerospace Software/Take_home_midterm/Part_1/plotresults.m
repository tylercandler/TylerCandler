function [] = plotresults(tspan,theta1_0, v_theta1_0, theta2_0, v_theta2_0, m1, m2, l1, l2, hand_position, hand_angle1, hand_angle2)
%function to plot the results of a double pendulum to the app
%set Initial conditions to a vector
y0 = [theta1_0, v_theta1_0, theta2_0, v_theta2_0]

%Call simulate pendulum to integrate initial conditions
[t, out] = simulate_pendulum(tspan, y0, l1, l2, m1, m2);



% Extract Angular Data from ode45 output
theta1 = out(:,1);
omega1 = out(:,2);
theta2 = out(:,3);
omega2 = out(:,4);


%Convert to cartesian coordinates, x and z 
for i = 1: length(t)
    x_1(i,1) = l1*sin(theta1(i));
    z_1(i,1) = -l1*cos(theta1(i));
    x_2(i,1) = l1*sin(theta1(i)) + l2*sin(theta2(i));
    z_2(i,1) = -l1*cos(theta1(i)) -l2*cos(theta2(i));
end



%Plot position angles to app
plot(hand_position,x_1,z_1)
hold (hand_position,'on')
plot(hand_position,x_2,z_2)
xlabel(hand_position,"X [m]")
ylabel(hand_position,"Y [m]")
title(hand_position,"Pendulum Visualization")
legend(hand_position,'Bob 1','Bob 2')

%Plot bob 2 angles to app
plot(hand_angle1,t,theta1)
hold (hand_angle1,'on')
xlabel(hand_angle1,"Time [s]")
ylabel(hand_angle1,"Angle [rads]")
title(hand_angle1,"Bob 1 Angle vs Time")

%Plot bob 1 angles to app
plot(hand_angle2,t,theta2)
hold (hand_angle2,'on')
xlabel(hand_angle2,"Time [s]")
ylabel(hand_angle2,"Angle [rads]")
title(hand_angle2,"Bob 2 Angle vs Time")

end