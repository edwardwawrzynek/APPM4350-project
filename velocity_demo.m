function velocity_demo
    set(0,'defaultTextInterpreter','latex');
    close all;
    % speed of light in free space
    c = 3.00e8;
    % waveguide mode cutoff
    omega_c = 2*pi*6.557e9;
    % frequencies to plot
    omega = omega_c:10e6:5*omega_c;

    vp = c ./ sqrt(1 - omega_c.^2 ./ omega.^2);
    vg = c.*sqrt(1 - omega_c.^2 ./ omega.^2);

    plot(omega ./ omega_c, vp, 'DisplayName', "v_p (phase velocity)");
    hold on;
    plot(omega ./ omega_c, vg, "DisplayName", "v_g (group velocity)");
    plot(omega ./ omega_c, c * ones(1, numel(omega)), 'k--', "HandleVisibility", "off");
    grid on;
    ylim([0 2*c]);
    xlabel("f/f_c");
    ylabel("velocity [m/s]");
    title("Phase and Group Velocity");
    legend();
end