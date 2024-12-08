clc; clear; close all;

%% Constants
mu0 = 4*pi*10^-7;
eps0 = 8.854*10^-12;
c0 = 1/sqrt(mu0*eps0);

%% Parameters
a = 22.86*10^-3; %0.9 inches
b = 10.16*10^-3; %0.4 inches
length = 14*10^-2; %5.5 inches


%% Extraction
filename = '14cm_line.s2p'; %input s2p file name
S = sparameters(filename);

s11 = rfparam(S, 1, 1);
s12 = rfparam(S, 1, 2);
s21 = rfparam(S, 2, 1);
s22 = rfparam(S, 2, 2); % Extract s-parameters from s2p file

freq = (S.Frequencies);
step = freq(2)-freq(1);

raw_phase = atan2(imag(s21),real(s21)); %phase from s2p
phase = unwrap(raw_phase); %unwrap phase
phase_deg = rad2deg(phase);

%% Math
m=1;
n=0; %define modes

k = (m*pi/a)^2+(n*pi/b)^2;
fc = (c0/2)*sqrt((m/a)^2+(n/b)^2);

calc_beta = (2*pi*freq./c0)*sqrt(1-(1/freq.^2)*((c0^2)*(k)/(4*pi^2)));
meas_beta = phase*-1/length; %calculate phase coefficient

v_phase_m = 2*pi*freq/meas_beta(:,1);
v_phase_c = 2*pi*freq/calc_beta(:,1); %calculate phase velocity

v_group_m = 1./(gradient(meas_beta(:))/step);
% v_group_m = sqrt(c0^2* meas_beta.^2+(2*pi*freq).^2)*(1/(meas_beta.*c0^2));
v_group_c = c0*sqrt(1-(fc./freq).^2); %calculate group velocity

%% Plotting
figure();
grid on
hold on

plot(freq/10^9,v_phase_m(:,201)/10^8,'LineWidth',1, color='#0072BD')
plot(freq/10^9,v_phase_c(:,201)/10^8,'LineWidth',1,color='#7E2F8E')
legend('Measured Phase Velocity','Calculated Phase Velocity','Interpreter', ...
    'latex','Location','southeast')

xlim([6.555 13])
ylim([0 4])
%ylim([-190 190])
title ("Phase Velocity of Waveguide",'Interpreter','latex')
xlabel('Frequency [GHz]','Interpreter','latex')
ylabel('Phase Velocity [$10^8$ m/s]','Interpreter','latex')
% ylabel('Angle [$^{\circ}]$','Interpreter','latex')

figure();
grid on
hold on
% 
plot(freq/10^9,v_group_m/10^8,'LineWidth',1, color='#0072BD')
plot(freq/10^9,v_group_c/10^8,'LineWidth',1,color='#7E2F8E')
legend('Measured Group Velocity','Calculated Phase Velocity','Interpreter', ...
    'latex','Location','southeast')

xlim([6.555 13])
ylim([0 4])
title ("Phase Velocity of Waveguide",'Interpreter','latex')
xlabel('Frequency [GHz]','Interpreter','latex')
ylabel('Group Velocity [$10^8$ m/s]','Interpreter','latex')
% ylabel('Angle [$^{\circ}]$','Interpreter','latex')
