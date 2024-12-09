% program to numerically calculate dispersion down waveguide
% speed of light in free space
c = 3.00e8;
% waveguide mode cutoff
omega_c = 2*pi*6.557e9;
% positions to sample at
x_samples = [0 0.3 0.9 1.6];
% pulse width
T = 1e-11;
% sampling period
samp_period = 1e-13;
omega_samp = 2*pi/samp_period;
% time sample over which to evaluate function
t = -1e-8:samp_period:1e-8;
% sample pulse waveform
u0 = cos(2*pi*10e9*t).*sin(2*pi*1e9*t).*(t > 0 & t < 0.5e-9);

% perform convolution numerically via fft
U0 = fftshift(fft(u0));
omega = omega_samp/numel(t)*(-numel(t)/2:numel(t)/2-1);
% apply dispersion relation
beta = sign(omega) .* sqrt((omega.^2 - omega_c.^2)./(c.^2));
% perform propogation
for x = x_samples
    U1 = U0 .* 0.5 .* (exp(-1i.*beta.*x) + exp(1i.*beta.*x)) .* (abs(omega) > omega_c);

    u1 = ifft(ifftshift(U1));
    plot(t*1e9, u1, "DisplayName", "x=" + x + " m");
    grid on;
    hold on;
    legend();
    xlabel("Time [ns]");
    ylabel("u");
end


