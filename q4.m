%%
V = 20; %volt
Kt = 0.46; % Nm/A
Kb = Kt; % Nm/A
R = 0.76; %ohm
L = 0;
J = 0.07; 
c = 0.5;
ca = 0.1;
cn = 0.075;
m = 0.8;
l = 0.3;

Wmax = 2*pi;
Wmax_dot = Wmax/0.2; 
Tmax = J*Wmax_dot+ca*(30/3.6)^2+c*Wmax+m*l*25;

% N=1, n=0.98
subplot(2,2,1)
N=1;n=0.98;
Wm_free = V/Kb;
WL_free = (1/N)*Wm_free;
WL = linspace(0,WL_free,1000);
TL = (Kt*n*N/(R))*(V-(Kt*N*WL))
plot(WL,TL,'LineWidth',2)
hold on
%plot(Wmax, Tmax, 'x')
title('N=1, n=0.98')
xlabel('w [rad/sec]')
ylabel('T [N*m]')
xlim([0, 45])
ylim([0, 80])
grid on

% N=3, n=0.95
subplot(2,2,2)
N=3;n=0.95;
Wm_free = V/Kb;
WL_free = (1/N)*Wm_free;
WL = linspace(0,WL_free,1000);
TL = (Kt*n*N/(R))*(V-(Kt*N*WL))
plot(WL,TL,'LineWidth',2)
hold on
%plot(Wmax, Tmax, 'x')
title('N=3, n=0.95')
xlabel('w [rad/sec]')
ylabel('T [N*m]')
xlim([0, 45])
ylim([0, 80])
grid on

% N=5, n=0.9
subplot(2,2,3)
N=5;n=0.9;
Wm_free = V/Kb;
WL_free = (1/N)*Wm_free;
WL = linspace(0,WL_free,1000);
TL = (Kt*n*N/(R))*(V-(Kt*N*WL))
plot(WL,TL,'LineWidth',2)
hold on
%plot(Wmax, Tmax, 'x')
title('N=5, n=0.9')
xlabel('w [rad/sec]')
ylabel('T [N*m]')
xlim([0, 45])
ylim([0, 80])
grid on

% N=7, n=0.88
subplot(2,2,4)
N=7;n=0.88;
Wm_free = V/Kb;
WL_free = (1/N)*Wm_free;
WL = linspace(0,WL_free,1000);
TL = (Kt*n*N/(R))*(V-(Kt*N*WL))
plot(WL,TL,'LineWidth',2)
hold on
%plot(Wmax, Tmax, 'x')
title('N=7, n=0.88')
xlabel('w [rad/sec]')
ylabel('T [N*m]')
xlim([0, 45])
ylim([0, 80])
grid on

%% working with N=3, n=0.95
% find a working point

N=3;n=0.95;

for V = [20 15 10 5]
    Wm_free = V/Kb;
    WL_free = (1/N)*Wm_free;
    WL = linspace(0,WL_free,1000);
    TL = (Kt/(R))*(V-(Kt*N*WL));
    plot(WL,TL,'LineWidth',2)
    hold on
end
WL = linspace(0,15,1000);
plot(WL,c*WL,'LineWidth',2)
title('N=3, n=0.95')
xlabel('w [rad/sec]')
ylabel('T [N*m]')
ylim([0, 20])
grid on

%% find open loop transfer function
%clear all
%syms R J c cn Kt N n Kb
s = tf('s');
V = 20; %volt
Kt = 0.46; % Nm/A
Kb = Kt; % Nm/A
R = 0.76; %ohm
L = 0;
J = 0.07; 
c = 0.5;
ca = 0.1;
cn = 0.075;
m = 0.8;
l = 0.3;
N=3;n=0.95;

P1 = 1/(R);
P2 = 1/(J*s^2+c*s+cn);

%Peq = ((P1*P2*Kt*N*n)/(1+Kb*N*P1*P2*Kt*N*n))*(1/s);
Pol = ((Kt*N*n)/(J*R))/(s^3+(R*c)/(J*R)*s^2+(1/(J*R)*(Kb*Kt*n*N^2)*s)) 
margin(Pol)

Kp = 2;
Ti = 3.4;
C_PI = Kp*(1+(1/(Ti*s)));
margin(C_PI*Pol)
%rlocus(C_PI*Pol)

t = out.theta.time;
theta = out.theta.signals.values;
plot(t, theta,'LineWidth', 2)
xlabel('time [sec]')
ylabel('theta [rad]')
%hold on
grid on

t = out.u.time;
u = out.u.signals.values;
plot(t, u,'LineWidth', 2)
xlabel('time [sec]')
ylabel('voltage [volt]')
%hold on
grid on

t = out.TL_sim.time;
TL_sim = out.TL_sim.signals.values;
%plot(t, TL_sim,'LineWidth', 2)
%xlabel('time [sec]')
%ylabel('TL [Nm]')
%hold on
grid on

t = out.WL_sim.time;
WL_sim = out.WL_sim.signals.values;
%plot(t, WL_sim,'LineWidth', 2)
%xlabel('time [sec]')
%ylabel('WL [rad/sec]')
%hold on
grid on

maxTL = max(TL_sim)
maxWL = max(WL_sim)

N=3;n=0.95;
Wm_free = V/Kb;
WL_free = (1/N)*Wm_free;
WL = linspace(0,WL_free,1000);
TL = (Kt*n*N/(R))*(V-(Kt*N*WL))
plot(WL,TL,'LineWidth',2)
hold on
plot(maxWL, maxTL,'o')
title('N=3, n=0.95')
xlabel('WL [rad/sec]')
ylabel('TL [N*m]')
xlim([0, 45])
ylim([0, 80])
grid on
hold off
