function dispersion_demo
    set(0,'defaultTextInterpreter','latex');
    close all;
    % speed of light in free space
    c = 3.00e8;
    % waveguide mode cutoff
    omega_c = 2*pi*6.557e9;
    % frequencies to plot
    omega1 = 2*pi*10e9;
    omega2 = 2*pi*10.5e9;

    % compute phase coefficients
    beta1 = sqrt((omega1^2 - omega_c^2)/(c^2));
    beta2 = sqrt((omega2^2 - omega_c^2)/(c^2));
    
    % compute parameters for
    omega = omega1;
    domega = omega2-omega1;

    beta = beta1;
    dbeta = beta2 - beta1;

    
    % time steps at which to plot
    times = [0 0.3e-9 0.6e-9];
    % spacial dimension over which to plot
    x = -0.4:0.001:0.4;
    
    plot_index = 1;
    for t = times
        % evaluate
        carrier = cos((beta+dbeta/2).*x - (omega + domega/2).*t);
        envelope = cos(dbeta/2.*x - domega/2.*t);
        u = 2*carrier.*envelope;

        carrier_vel = (2*omega + domega) / (2*beta + dbeta);
        envelope_vel = (domega/dbeta);
        
        % plot
        ax = subplot(numel(times),1,plot_index);
        plot(x, u, 'b', "DisplayName", "u");
        hold on;
        plot(x, 2*envelope, 'k--', "DisplayName", "envelope");
        plot(x, -2*envelope, 'k--', "HandleVisibility", "off");

        plot([carrier_vel*t], [2*cos(dbeta/2*carrier_vel*t - domega/2*t)], 'r.', 'MarkerSize', 15, "DisplayName", "fixed phase");
        
        grid on;

        ylim([-2.2 2.2]);
        xlim([min(x) max(x)]);

        if plot_index == 1
            legend();
            xlabel(["x[m]", "(a) t = 0ns"]);
        elseif plot_index == 2
            xlabel(["x[m]", "(b) t = 0.3ns"]);
        else
            xlabel(["x[m]", "(c) t = 0.6ns"]);
        end
        ylabel("u");

        plot_index = plot_index + 1;
    end
end