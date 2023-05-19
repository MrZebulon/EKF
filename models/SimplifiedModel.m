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
        
        % Time has been taken into account within the biases/noises/drifts
        % units = m/s
        accel_bias = [0.0029394893123127377, -0.0009818539510592773, 0.0028762637247315066];
        accel_noise = 4.7358607479753485e-09;
        accel_drift = 3.314312818032142e-10;

        % units = m
        baro_bias = 399.23657624056926;
        baro_noise = 0.014769875002697693;
        baro_drift = 6.282361435771973e-05;
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
                0, 0, 0,    0, 0, 0,         -1, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, -1, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, -1,      0
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
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0
                0, 0, 0,    dvxCov, 0, 0,      0, 0, 0,      0
                0, 0, 0,    0, dvyCov, 0,      0, 0, 0,      0
                0, 0, 0,    0, 0, dvzCov,      0, 0, 0,      0
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0
                0, 0, 0,    0, 0, 0,           0, 0, 0,      0];
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
            vel_delta_bias_sigma = scale_var.* obj.accel_drift;
            pos_delta_bias_sigma = scale_var.* obj.baro_drift;

            Qs = diag([obj.additive_noise.*ones(1,6), vel_delta_bias_sigma*ones(1,3), pos_delta_bias_sigma*ones(1,1)]);
            w = scale_var.*[obj.accel_noise*ones(1,3)];
        end

        function R = get_R_matrix(obj)
            R = obj.baro_noise;
        end
    end
end

