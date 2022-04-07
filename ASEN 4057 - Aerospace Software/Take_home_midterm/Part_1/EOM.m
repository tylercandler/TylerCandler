%% Equations of Motion function handle
function dydt = EOM(t,g,L1,L2,m1,m2,state) %state = [theta1 v_theta1, theta2, v_theta2]'
%Extract data from system state
theta1 = state(1);
v_theta1 = state(2);
theta2 = state(3);
v_theta2 = state(4);

%Substitutions as per documentation
a = (m1+m2)/L1;
b = m2*L2*cos(theta1-theta2);
c = m2*L2;
d = m2*L2*cos(theta1-theta2);
e = m2*L2*v_theta2^2*sin(theta1-theta2) + g*(m1+m2)*sin(theta1);
f = m2*L1*v_theta1^2*sin(theta1-theta2) + m2*g*sin(theta2);

%Theta1_ dotdot and Theta2_dotdot
a_theta2 = -(f + (d*e)/a) / (c - (d*b)-a);

a_theta1 = (-e - (b * a_theta2))/a;

%Output dydt
dydt = [v_theta1, a_theta1, v_theta2, a_theta2]';

end