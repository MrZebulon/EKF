classdef SimplifiedModel < BaseModel
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
    properties (Constant)
        additive_noise = 1e-8;
        
        % units = g (9.81 m/s)
        accel_bias = [-0.2668811881 0.3305544554 0.07695049505];
        accel_noise = 7.05E-04;
        accel_bias_noise =  6.89e-4;

        % units = hPa
        baro_bias = 970.5698218;
        baro_noise = 1.52e-5;
        baro_bias_noise = 2.98e-7;

        baro_measurement_uncertainty = 0.1;
    end

    methods
        function [x_init, P_init] = get_init_state(obj)
            x_init = [
                0
                0
                0
                0
                0
                0
                obj.accel_bias(1)
                obj.accel_bias(2)
                obj.accel_bias(3)
                obj.baro_bias];

            P_init = diag(1e-9 * ones(1, 10));
        end

        function dx = get_delta_x(obj, x, u)
            dt = obj.dt;
            dx = [
                x(4) * dt
                x(5) * dt
                x(6) * dt
                (u(1) + x(7)) * dt
                (u(2) + x(8)) * dt
                (u(3) + x(9)) * dt
                0
                0
                0
                0];
        end

        function F = get_F_matrix(obj, x, u)
            F = [...
                0, 0, 0,    1, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 1, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 1,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    1, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 1, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 1,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0
                0, 0, 0,    0, 0, 0,    0, 0, 0,    0];
        end

        function G = get_G_matrix(obj, x, u, w)
            G = [...
                0, 0, 0     0,   0, 0, 0,   0
                0, 0, 0     0,   0, 0, 0,   0
                0, 0, 0     0,   0, 0, 0,   0
                1, 0, 0     0,   0, 0, 0,   0
                0, 1, 0     0,   0, 0, 0,   0
                0, 0, 1     1,   0, 0, 0,   0
                0, 0, 0     0,   1, 0, 0,   0
                0, 0, 0     0,   0, 1, 0,   0
                0, 0, 0     0,   0, 0, 1,   0
                0, 0, 0     0,   0, 0, 0,   1];
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

