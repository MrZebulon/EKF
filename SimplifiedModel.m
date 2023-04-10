classdef SimplifiedModel < Model
    %
    % State vector definition
    %   1-3 : position
    %   4-6 : speed
    %   7-9 : accel bias
    %   10 : baro bias
    %
    % Control vector definition
    %   1-3 : acceleration
    %
    % Noise vector definition
    %   1-3 : acceleration noise
    %   4 : barometer noise
    %   5-7 : acceleration bias noise
    %   8 : barometer bias noise
    properties
        additive_noise = 1e-8;

        accel_noise = 7.05E-04;
        accel_bias_noise =  6.89e-4;

        baro_noise = 1.52e-5; % units = hPa
        baro_bias_noise = 2.98e-7; % units = hPa
        baro_measurement_uncertainty = 0.1;
    end

    methods
        function x_new = get_next_state(obj, x, u)
            dt = obj.dt;
            x_new = [
                x(4) * dt
                x(5) * dt
                x(6) * dt + x(10)
                u(1) * dt + x(7)
                u(2) * dt + x(8)
                u(3) * dt + x(9)
                x(7)
                x(8)
                x(9)
                x(10)];
        end

        function F = get_F_matrix(obj, x, u)
            F = [...
                0, 0, 0,    1, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 1, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 1,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0];
        end

        function G = get_G_matrix(obj, x, u, w)
            cov_acc_n = w(1);
            cov_acc_e = w(2);
            cov_acc_d = w(3);
            cov_bias_acc_n = w(4);
            cov_bias_acc_e = w(5);
            cov_bias_acc_d = w(6);
            cov_baro = w(7);
            cov_bias_baro = w(8);

            G = [...
                0, 0, 0             0, 0, 0,                0, 0
                0, 0, 0             0, 0, 0,                0, 0
                0, 0, 0             0, 0, 0,                cov_baro, 0
                cov_acc_n, 0, 0     0, 0, 0,                0, 0
                0, cov_acc_e, 0     0, 0, 0,                0, 0
                0, 0, cov_acc_d     cov_bias_acc_n, 0, 0,   0, 0
                0, 0, 0             0, cov_bias_acc_e, 0,   0, 0
                0, 0, 0             0, 0, cov_bias_acc_d,   0, 0
                0, 0, 0             0, 0, 0,                0, 0
                0, 0, 0             0, 0, 0,                0, cov_bias_baro];
        end

        function z_hat  = get_measurement_estimate(obj, x)
            z_hat = x(3);
        end

        function H = get_H_matrix(obj)
            H = [...
                0, 0, 1,    0, 0, 0,    0, 0, 0,    0];
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            scale_var = 0.5*(1./(Fs.^2));
            vel_delta_bias_sigma = scale_var.* obj.accel_bias_noise;
            pos_delta_bias_sigma = scale_var.* obj.baro_bias_noise;

            Qs = diag([obj.additive_noise.*ones(1,3), vel_delta_bias_sigma*ones(1,3), obj.additive_noise.*ones(1,1), pos_delta_bias_sigma*ones(1,1)]);
            w = scale_var.*[obj.accel_noise*ones(1,3), obj.accel_bias_noise*ones(1,3), obj.baro_noise*ones(1,1), obj.baro_bias_noise*ones(1, 1)];
        end

        function R = get_R_matrix(obj)
            R = obj.baro_measurement_uncertainty;
        end
    end
end

