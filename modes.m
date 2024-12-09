clc; clear; close all;

%% Constants
mu0 = 4*pi*10^-7;
eps0 = 8.854*10^-12;
c0 = 1/sqrt(mu0*eps0);

%% Parameters
a = 22.86*10^-3; %0.9 inches
b = 10.16*10^-3; %0.4 inches
length = 14*10^-2; %5.5 inches
freq = 8.5*10^9; %8.5GHz, well within band
w = freq*2*pi;
er = 1; %permittivity of inside material
k = w/(c0*sqrt(er));

res = 25; %resolution of image


%% TE Mode
m = 2;
n = 0;
A = 1; % magnitude
eigenm = m*pi/a;
eigenn = n*pi/b;

x = linspace(0,a,res);
test = x;
y = linspace(0,b,res);
z = 0;

[x,y] = meshgrid(x,y);

kc = sqrt(eigenm^2 + eigenn^2);
beta = sqrt(k^2 - kc^2);

%Electric Field Equations
Ex = 1j*w*mu0*n*pi/(b*kc^2)*A*cos(eigenm.*x).*sin(eigenn.*y).*exp(-1j*beta*z);
Ey= -1j*w*mu0*m*pi/(a*kc^2)*A*sin(eigenm.*x).*cos(eigenn.*y).*exp(-1j*beta*z);

%Magnetic Field Equations
Hx = 1j*beta*m*pi/(a*kc^2)*A*sin(eigenm.*x).*cos(eigenn.*y).*exp(-1j*beta*z);
Hy = 1j*beta*n*pi/(b*k^2)*A*cos(eigenm.*x).*sin(eigenn.*y).*exp(-1j*beta*z);

%% TE Plotting
% Front View
figure();
hold on
quiver(x.*10^3,y.*10^3,imag(Ex),imag(Ey),1.1,'color','#A2142F')
xlim([0 a*10^3])
ylim([0 b*10^3])
%title ("E-Field $TM_{"+string(m)+string(n)+"}$"','Interpreter','latex')
% xlabel('x-position [mm]','Interpreter','latex')
% ylabel('y-position [mm]','Interpreter','latex')

quiver(x.*10^3,y.*10^3,real(Hx),real(Hy),'color','#0072BD')

title ("$TE_{"+string(m)+string(n)+"}$",'Interpreter','latex')
xlabel('x-position [mm]','Interpreter','latex')
ylabel('y-position [mm]','Interpreter','latex')

legend('E-Field','H-Field')

%% TM Mode
m = 2;
n = 1;
B = 1; % magnitude
eigenm = m*pi/a;
eigenn = n*pi/b;

x = linspace(0,a,res);
y = linspace(0,b,res);
z=0;
[x,y] = meshgrid(x,y);


kc = sqrt(eigenm^2 + eigenn^2);
beta = sqrt(k^2 - kc^2);

Ex = -beta*m*pi*(a*kc^2)*B*cos(eigenm*x).*sin(eigenn.*y);
Ey= -beta*m*pi*(b*kc^2)*B*sin(eigenm*x).*cos(eigenn.*y);
Ez = B*sin(eigenm.*x).*sin(eigenn*y).*exp(-1j*beta*z);

Hx = w*er*n*pi/(b*kc^2)*B*sin(eigenm*x).*cos(eigenn.*y);
Hy = -w*er*n*pi/(a*k^2)*B*cos(eigenm*x).*sin(eigenn.*y);
Hz = 0;

%% TM Plotting
% Front View
figure();
hold on
quiver(x.*10^3,y.*10^3,imag(Ex),imag(Ey),'color','#A2142F')
% xlim([0 a*10^3])
% ylim([0 b*10^3])
% title ("E-Field $TM_{"+string(m)+string(n)+"}$",'Interpreter','latex')
% xlabel('x-position [mm]','Interpreter','latex')
% ylabel('y-position [mm]','Interpreter','latex')

quiver(x.*10^3,y.*10^3,real(Hx),real(Hy),'color','#0072BD')
xlim([0 a*10^3])
ylim([0 b*10^3])
title ("$TM_{"+string(m)+string(n)+"}$",'Interpreter','latex')
xlabel('x-position [mm]','Interpreter','latex')
ylabel('y-position [mm]','Interpreter','latex')
legend('E-Field','H-Field')




