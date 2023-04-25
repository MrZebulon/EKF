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
        
        % units = m/s
        accel_bias = [-2.684265624256846 3.3420780321046384 9.018026218787151];
        accel_noise = 2;
        accel_bias_noise = 0.0003718319715151406;

        % units = m
        baro_bias = 361.3487972164834;
        baro_noise = 0.0011529662809109404;
        baro_bias_noise = 2.8486463440220755e-06;

        baro_measurement_uncertainty = 1;
    end

    methods
        function [x_init, P_init] = get_init_state(obj)
            x_init = [
                0
                0
                obj.baro_bias

                0
                0
                0
                
                obj.accel_bias(1)
                obj.accel_bias(2)
                obj.accel_bias(3)

                obj.baro_bias];

            P_init = diag(1e-9 * ones(1, 10));
        end

        function x_new = compute_x_new(obj, x, u)
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

            x_new = x + dx;
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

