function output_derivative = ODEFUN(t,input)
%% Givens & Variable explanations
G=6.67*10^(-11);
mM = 7.34767309*10^22; % mass of moon in kg
mE = 5.97219*10^24; % mass of Earth kg
mS = 28833; % mass of spacecraft in kg
rM = 1737100; %radius of moon in m
rE = 6371000; %r_adius of the Earth in m
% d_EM=384403000; %distance from Earth to moon in m

%Extract from input vec
xS = input(1);
yS = input(2);
vSx = input(3);
vSy = input(4);
xM = input(5);
yM = input(6);
vMx = input(7);
vMy = input(8);
xE = input(9);
yE = input(10);
vEx = input(11);
vEy = input(12);



%Force and Distance relations
%Force from Earth to Moon (F_EM)
d_EM = sqrt((xM-xE)^2 + (yM-yE)^2);
% (xM-xE)^2
% (xM-xE)^2

F_EMx=(G*mE*mM*(xE-xM))/(d_EM^3); % force on Moon by Earth in the x-direction
F_EMy=(G*mE*mM*(yE-yM))/(d_EM^3); % force on Moon by Earth in the y-direction

%Force from Satellite to Moon (F_SM)
d_SM = sqrt((xS-xM)^2 + (yS-yM)^2);

F_MSx=(G*mS*mM*(xM-xS))/(d_SM^3); % force on Satellite by Moon in the x-direction
F_MSy=(G*mS*mM*(yM-yS))/(d_SM^3); % force on Satellite by Moon in the y-direction


%Force from Satellite to Earth (F_ES)
d_ES = sqrt((xS-xE)^2 + (yS-yE)^2);

F_ESx=(G*mS*mE*(xE-xS))/(d_ES^3); % force on Satellite by Moon in the x-direction
F_ESy=(G*mS*mE*(yE-yS))/(d_ES^3); % force on Satellite by Moon in the y-direction



% dAB=sqrt((xA-xB)^2+(yA-yB)^2); %distance between bodies A and B
% FAB=(G*mA*mB)/(dAB^2); % general force on B by A
% FABx=(G*mA*mB*(xA-xB))/(dAB^3); % force on B by A in the x-direction
% FABy=(G*mA*mB*(yA-yB))/(d_AB^3); % force on B by A in the y-direction
% FBAx=(G*m_a*m_b*(xB-xA))/(d_AB^3); % force on A by B in the x-direction
% FBAy=(G*m_a*m_b*(yB-yA))/(d_AB^3); % force on A by B in the y-direction


%Newton's second law and acceleration
aSx=(F_MSx+F_ESx)/mS;
aSy=(F_MSy+F_ESy)/mS;
aMx=(F_EMx-F_MSx)/mM;
aMy=(F_EMy-F_MSy)/mM;
aEx=0;
aEy=0;
% aEx = (F_EMx + F_ESx)/mE;
% aEy = (F_EMy + F_ESy)/mE;

 %derivative acceleration relations
 %calculated via relative law of gravitation

%  aSx=dvSx/dt;
%  aSy=dvSy/dt;
%  aMx=dvMx/dt;
%  aMy=fvMy/dt;
%  aEx=dvEx/dt;
%  aEy=dvEy/dt;

 %derivative velocity relations 
 %directly from input arguments
% vSx=dxS/dt;
% vSy=dyS/dt;
% vMx=dxM/dt;
% vMy=dyM/dt;
% vEx=dxE/dt;
% vEy=dyE/dt;

% deltaT
% if length(t)==1;
%     deltaT=0;
% else 
%     deltaT=t(end)-t(end-1);
% end

output_derivative = [vSx, vSy, aSx, aSy, vMx, vMy, aMx, aMy, vEx, vEy, aEx, aEy]';
% 



end