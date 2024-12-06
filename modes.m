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
m = 1;
n = 1;
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

Ex = 1j*w*mu0*n*pi/(b*kc^2)*A*cos(eigenm.*x).*sin(eigenn.*y).*exp(-1j*beta*z);
Ey= -1j*w*mu0*m*pi/(a*kc^2)*A*sin(eigenm.*x).*cos(eigenn.*y).*exp(-1j*beta*z);
% Ez = 0;

Hx = 1j*beta*m*pi/(a*kc^2)*A*sin(eigenm.*x).*cos(eigenn.*y).*exp(-1j*beta*z);
Hy = 1j*beta*n*pi/(b*k^2)*A*cos(eigenm.*x).*sin(eigenn.*y).*exp(-1j*beta*z);
% Hz = A*cos(eigenm.*x).*cos(eigenn.*y).*exp(-1j*beta.*z);
%% TE Plotting
% Front View
figure();
quiver(x.*10^3,y.*10^3,imag(Ex),imag(Ey))
xlim([0 a*10^3])
ylim([0 b*10^3])
title ('E-Field $TE_{10}$','Interpreter','latex')
xlabel('x-position [mm]','Interpreter','latex')
ylabel('y-position [mm]','Interpreter','latex')

figure();
quiver(x.*10^3,y.*10^3,imag(Hx),imag(Hy))
xlim([0 a*10^3])
ylim([0 b*10^3])
title ('H-Field $TE_{10}$','Interpreter','latex')
xlabel('x-position [mm]','Interpreter','latex')
ylabel('y-position [mm]','Interpreter','latex')

%% TM Mode
m = 2;
n = 0;
B = 1; % magnitude
eigenm = m*pi/a;
eigenn = n*pi/b;

x = linspace(0,a,res);
y = linspace(0,b,res);
z = linspace(0,length,res);

[x,y] = meshgrid(x,y);


kc = sqrt(eigenm^2 + eigenn^2);
beta = sqrt(k^2 - kc^2);

Ex = -beta*m*pi*(a*kc^2)*B*cos(eigenm*x).*sin(eigenn.*y);
Ey= -beta*m*pi*(b*kc^2)*B*sin(eigenm*x).*cos(eigenn.*y);
Ez = B*sin(eigenm.*x).*sin(eigenn*y).*exp(-1j*beta*z);

Hx = w*er*n*pi/(b*kc^2)*B*sin(eigenm*x).*cos(eigenn.*y).*exp(-1j*beta*z);
Hy = -w*er*n*pi/(a*k^2)*B*cos(eigenm*x).*sin(eigenn.*y).*exp(-1j*beta*z);
Hz = 0;

% TM Plotting
% Front View
figure();
quiver(x.*10^3,y.*10^3,imag(Ex),imag(Ey))
xlim([0 a*10^3])
ylim([0 b*10^3])
title ('E-Field $TM_{10}$','Interpreter','latex')
xlabel('x-position [mm]','Interpreter','latex')
ylabel('y-position [mm]','Interpreter','latex')

figure();
quiver(x.*10^3,y.*10^3,real(Hx),real(Hy))
xlim([0 a*10^3])
ylim([0 b*10^3])
title ('H-Field $TM_{10}$','Interpreter','latex')
xlabel('x-position [mm]','Interpreter','latex')
ylabel('y-position [mm]','Interpreter','latex')



