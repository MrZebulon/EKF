classdef TranslationRotationModelUKF < BaseUKFModel
    %
    % State vector definition
    %   1-3 : position
    %   4-6 : velocity
    %   7-10 : orientation
    %   11-13 : accel bias
    %   14-16 : gyro bias
    %   17-19 : earth permanent field in inertial frame
    %   20-22: magneto_bias
    %
    % Control vector definition
    %   1-3 : acceleration
    %   4-6 : angular velocity
    %
    % Noise vector definition
    %   1-3 : accel noise
    %   4-7 : 0, gyro noise
    %
    % Obs vector definition
    %   1 : baro
    %   2-4 : magneto

    properties (Constant)
        Nx = 22;
        Nu = 6;
        Nw = 7;
    end
    properties (Constant)
        additive_noise = 1e-8;

        % Time has been taken into account within the biases/noises/drifts
        % units = m/s
        accel_bias = [-0.00022021696252465473, 1.972386587771194e-07, -0.0032535502958579853];
        accel_noise = 4.7358607479753485e-09;
        accel_drift = 3.314312818032142e-10;

        %units = rad/s
        gyro_bias = [4.9309664694280084e-05, -1.9723865877712034e-05, 0];
        gyro_noise = 1.0102261028815603e-08;
        gyro_drift = 3.9979150848914756e-10;

        % units = m
        baro_bias = 402.42156049763133;
        baro_noise = 0.014769875002697693;

        % units = uT
        magneto_bias = [100, 100, 100];
        magneto_noise = 0.09;
        magneto_drift =  1e-10;
        magneto_earth_field_drift =  1e-15;
        magneto_earth_field = [20.31736686390534, -6.346311637080932, -66.59784023668604]

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

                1
                0
                0
                0

                obj.accel_bias(1)
                obj.accel_bias(2)
                obj.accel_bias(3)

                obj.gyro_bias(1)
                obj.gyro_bias(2)
                obj.gyro_bias(3)

                obj.magneto_earth_field(1)
                obj.magneto_earth_field(2)
                obj.magneto_earth_field(3)

                obj.magneto_bias(1)
                obj.magneto_bias(2)
                obj.magneto_bias(3)];

            P_init = diag(1e-9 * ones(1, obj.Nx));
        end

        function x_new = state_prediction_function(obj, x, u)
            dt = obj.dt;

            x_new = [
                x(1:3) + x(4:6) * dt;
                x(4:6) + Utils.rotmat(x(7:10)') * (dt*u(1:3)' - x(11:13));
                Utils.normalize(Utils.mult_quat(x(7:10),[1.0; (dt.*u(4:6)'-x(14:16))/2]));
                zeros(3);
                zeros(3);
                zeros(3);
                zeros(3);
            ];

        end

        function Q = get_Q_matrix(obj, x, u, w)

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            G = zeros(obj.Nx, obj.Nw);

            G(4:6, 1:3) = quat2rotm([q0, q1, q2, q3]);
            G(7:10, 4:7) = 0.5 .*Utils.hamilton_product_as_matrix([q0, q1, q2, q3]);
            Q = G*diag([w(1:3), w(4:7)])*(G.');

        end

       function h = baro_obs_function(obj, x)
            h = x(3) + obj.baro_bias;
        end

        function R = baro_R_matrix(obj, x)
            R = obj.baro_noise;
        end

        function h = mag_obs_function(obj, x)
            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            magNavX = x(17); magNavY = x(18); magNavZ = x(19);
            magBiasX = x(20); magBiasY = x(21); magBiasZ = x(22);

            mx = magBiasX + magNavX*(q0^2 + q1^2 - q2^2 - q3^2) - magNavZ*(2*q0*q2 - 2*q1*q3) + magNavY*(2*q0*q3 + 2*q1*q2);
            my = magBiasY + magNavY*(q0^2 - q1^2 + q2^2 - q3^2) + magNavZ*(2*q0*q1 + 2*q2*q3) - magNavX*(2*q0*q3 - 2*q1*q2);
            mz = magBiasZ + magNavZ*(q0^2 - q1^2 - q2^2 + q3^2) - magNavY*(2*q0*q1 - 2*q2*q3) + magNavX*(2*q0*q2 + 2*q1*q3);
            
            h = [mx, my, mz];
        end

        function R = mag_R_matrix(obj, x)
            R = diag(obj.magneto_noise.*ones(1,3));
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            scale_var = 0.5*(1./(Fs.^2));
            accel_drift_sigma = scale_var.* obj.accel_drift;
            gyro_drift_sigma = scale_var.* obj.gyro_drift;

            Qs = diag([obj.additive_noise.*ones(1,10), accel_drift_sigma.*ones(1,3), gyro_drift_sigma.*ones(1,3), obj.magneto_earth_field_drift.*ones(1,3), obj.magneto_drift.*ones(1,3)]);
            w = scale_var.*[obj.accel_noise.*ones(1,3), 0, obj.gyro_noise.*ones(1,3)];
        end

    end
end
