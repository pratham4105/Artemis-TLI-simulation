# Artemis-TLI-simulation
% initial Data
RE = 6371e3
RM = 1737e3
dEM = 384400e3 %distance from center of earth to center of moon
g0E = 9.807 % m/s^2
g0M = 1.62 %m/s^2
muE = g0E*RE^2
TM = 2*pi*sqrt(dEM^3/muE) %seconds
TMdays = TM/60/60/24 % period in days
omM = 2*pi/TM % rad/s, angular speed of earth-moon rotating frame

% TLI orbit data
hp = 185e3 % m, perigee altitude (from surface of earth)
dam = 6500e3 % distance from moon surface to apogee
rp = RE + hp % perigee from center of earth
ra = dEM + RM + dam % apogee from center of the earth
a = (rp+ra)/2 % semi-major axis 
c = a - rp % location of foci
b = sqrt(a^2 - c^2) % semi-minor axis
vp = sqrt(muE*(2/rp-1/a)) % m/s, initial speed at perigee, vis-visa
va = sqrt(muE*(2/ra-1/a)) % m/s, speed at apogee, vis-visa
T = 2*pi*sqrt(a^3/muE) % s, orbital period of artemis
Thrs = T/60/60 % period of artemis in hours

thetaM_0 = 2*pi*((T/2)/TM) % radians, initial moon angle from horizon at t = 0

%initial condition vector
ic = [-rp;0;0;-vp] % [x;y;vx;vy]

%simulation time
time = [0 T]

%solving ODE
options = odeset('RelTol',1e-9);
[t,y] = ode45(@(t,y) inertialEOM(t,y,muE), time, ic, options)

% extracting and calculating values
x = y(:,1);
y_pos = y(:,2)
vx = y(:,3)
vy = y(:,4)

thetaMoon = omM*t - thetaM_0
xMoon = dEM*cos(thetaMoon)
yMoon = dEM*sin(thetaMoon)
dAM = sqrt((x-xMoon).^2 + (y_pos-yMoon).^2) % distance between artemis and moon

% figure 1 : artemis orbit
figure
plot(x/1000, y_pos/1000,'r','LineWidth',1.5)
hold on
theta = linspace(0,2*pi,200)
xEarth = RE*cos(theta)
yEarth = RE*sin(theta)
plot(xEarth/1000, yEarth/1000,'b', 'LineWidth', 1.5)
axis equal
grid on
xlabel('x (km)')
ylabel('y (km)')
title('Artemis TLI orbit - intertial frame')
legend('Artemis', 'Earth', 'Location', 'best')

% artemis moon distance figure
figure
plot(t/3600, dAM/1000, 'LineWidth', 1.5)
grid on
xlabel('Time (hours)')
ylabel('distance from artemis to moon (km)')
title('distance between artemis and moon')

% Equations of motion
function dydt = inertialEOM(t,y,muE)
    x = y(1)
    y_pos = y(2)
    vx = y(3)
    vy = y(4)

    r = sqrt(x^2 + y_pos^2);

    ax = -muE*x/r^3
    ay = -muE*y_pos/r^3

    dydt = [vx;vy;ax;ay];
end
