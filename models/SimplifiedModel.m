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
        accel_bias = [-0.026842656242568 0.033420780321046 -0.007947030636161];
        accel_noise = 1;
        accel_bias_noise = 2e-4;

        % units = m
        baro_bias = 361.3487972164834;
        baro_noise = 0.0011529662809109404;
        baro_bias_noise = 2.8486463440220755e-06;

        baro_measurement_uncertainty = 0.1;
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

                u(1)*dt - x(7)
                u(2)*dt - x(8)
                u(3)*dt - x(9)

                0
                0
                0
                
                0];

            x_new = x + dx;
        end
        
        function F = get_F_matrix(obj, x, u)
            % The F matrix is the Jacobian matrix of the state update
            % function f(x, u, w) with respect to the state vector x.
            F = [...
                0, 0, 0,    obj.dt, 0, 0,    0, 0, 0,      0
                0, 0, 0,    0, obj.dt, 0,    0, 0, 0,      0
                0, 0, 0,    0, 0, obj.dt,    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         1, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 1, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 1,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0];
        end

        function Q = get_Q_matrix(obj, x, u, w)
           
            dvxCov = w(1); 
            dvyCov = w(2);
            dvzCov = w(3); 
           
           Q = [...
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    dvxCov*obj.dt^2, 0, 0,      0, 0, 0,      0
                0, 0, 0,    0, dvyCov*obj.dt^2, 0,      0, 0, 0,      0
                0, 0, 0,    0, 0, dvzCov*obj.dt^2,      0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0];
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

            Qs = diag([obj.additive_noise.*ones(1,6), vel_delta_bias_sigma*ones(1,3), pos_delta_bias_sigma*ones(1,1)]);
            w = scale_var.*[obj.accel_noise*ones(1,3), obj.baro_noise*ones(1,1)];
        end

        function R = get_R_matrix(obj)
            R = obj.baro_measurement_uncertainty;
        end
    end
end

