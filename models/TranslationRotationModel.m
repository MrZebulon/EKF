classdef TranslationRotationModel < BaseModel
    %
    % State vector definition
    %   1-3 : position
    %   4-6 : velocity
    %   7-10 : orientation
    %   11-13 : accel bias
    %   14-16 : gyro bias
    %   17 : baro bias
    %
    % Control vector definition
    %   1-3 : acceleration
    %   4-6 : angular velocity 

    % Noise vector definition
    %   1-3 : acceleration noise
    %   4-6 : gyro noise
    %   7 : barometer noise
    %   8-10 : acceleration bias noise
    %   11-13 : gyro bias noise
    %   14 : barometer bias noise

    properties (Constant)
        additive_noise = 1e-8;
        
        % units = m/s
        accel_bias = [-2.6845463600451454 3.335545323081266 9.019774604966166];
        accel_noise = 0.0022501515160489754;
        accel_bias_noise = 0.0003718319715151406;

        %units = rad/s
        gyro_bias = [0.00583044566863878 -0.0005276800977633389 0.005892403831812206];
        gyro_noise = 0.0002976371540218863;
        gyro_bias_noise = 1.2221470219227153e-05;

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

                0
                1/sqrt(3)
                1/sqrt(3)
                1/sqrt(3)
                
                obj.accel_bias(1)
                obj.accel_bias(2)
                obj.accel_bias(3)

                obj.gyro_bias(1)
                obj.gyro_bias(2)
                obj.gyro_bias(3)

                obj.baro_bias];

            P_init = diag(1e-9 * ones(1, 17));
        end
        
        function x_new = compute_x_new(obj, x, u)
            dt = obj.dt;
            delta_acc = Utils.body_to_inertial_frame(x(7:10)', (u(1:3)' + x(14:16)));
            
            dx = [
                x(4) * dt
                x(5) * dt
                x(6) * dt

                delta_acc(1) * dt
                delta_acc(2) * dt
                delta_acc(3) * dt

                0 % delta_q must be computed separately
                0
                0
                0
                
                0
                0
                0

                0
                0
                0
                
                0];
            
            x_new = x + dx;
            x_new(7:10) = compact(normalize(quaternion(x(7:10)') * quaternion(u(4:6) + x(14:16)', 'rotvec')));
        end

        function F = get_F_matrix(obj, x, u)
            
            F = zeros(17, 17);

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            a_x = u(1); a_y = u(2); a_z = u(3);
            w_x = u(4); w_y = u(5); w_z = u(6);
            b_acc_x = x(11); b_acc_y = x(12); b_acc_z = x(13);
            b_gyro_x = x(14); b_gyro_y = x(15); b_gyro_z = x(16); 

            % dp^dot
            F(1:3, 4:6) = eye(3);
    
            %dv^dot
            F(4, :) = [0, 0, 0, 0, 0, 0, 2*conj(q2)*(a_z + b_acc_z) - 2*conj(q3)*(a_y + b_acc_y),                             2*conj(q2)*(a_y + b_acc_y) + 2*conj(q3)*(a_z + b_acc_z), 2*conj(q1)*(a_y + b_acc_y) - 4*conj(q2)*(a_x + b_acc_x) + 2*conj(q0)*(a_z + b_acc_z), 2*conj(q1)*(a_z + b_acc_z) - 2*conj(q0)*(a_y + b_acc_y) - 4*conj(q3)*(a_x + b_acc_x),           1 - 2*conj(q3)^2 - 2*conj(q2)^2, 2*conj(q1)*conj(q2) - 2*conj(q0)*conj(q3), 2*conj(q0)*conj(q2) + 2*conj(q1)*conj(q3), 0, 0, 0, 0];
            F(5, :) = [0, 0, 0, 0, 0, 0, 2*conj(q3)*(a_x + b_acc_x) - 2*conj(q1)*(a_z + b_acc_z), 2*conj(q2)*(a_x + b_acc_x) - 4*conj(q1)*(a_y + b_acc_y) - 2*conj(q0)*(a_z + b_acc_z),                             2*conj(q1)*(a_x + b_acc_x) + 2*conj(q3)*(a_z + b_acc_z), 2*conj(q0)*(a_x + b_acc_x) - 4*conj(q3)*(a_y + b_acc_y) + 2*conj(q2)*(a_z + b_acc_z), 2*conj(q0)*conj(q3) + 2*conj(q1)*conj(q2),           1 - 2*conj(q3)^2 - 2*conj(q1)^2, 2*conj(q2)*conj(q3) - 2*conj(q0)*conj(q1), 0, 0, 0, 0];
            F(6, :) = [0, 0, 0, 0, 0, 0, 2*conj(q1)*(a_y + b_acc_y) - 2*conj(q2)*(a_x + b_acc_x), 2*conj(q3)*(a_x + b_acc_x) + 2*conj(q0)*(a_y + b_acc_y) - 4*conj(q1)*(a_z + b_acc_z), 2*conj(q3)*(a_y + b_acc_y) - 2*conj(q0)*(a_x + b_acc_x) - 4*conj(q2)*(a_z + b_acc_z),                             2*conj(q1)*(a_x + b_acc_x) + 2*conj(q2)*(a_y + b_acc_y), 2*conj(q1)*conj(q3) - 2*conj(q0)*conj(q2), 2*conj(q0)*conj(q1) + 2*conj(q2)*conj(q3),           1 - 2*conj(q2)^2 - 2*conj(q1)^2, 0, 0, 0, 0];
            %dq^dot
            F(7, :) = [0, 0, 0, 0, 0, 0,               1/2, - b_gyro_x/2 - w_x/2, - b_gyro_y/2 - w_y/2, - b_gyro_z/2 - w_z/2, 0, 0, 0, -q1/2, -q2/2, -q3/2, 0];
            F(8, :) = [0, 0, 0, 0, 0, 0, b_gyro_x/2 + w_x/2,                 1/2,   b_gyro_z/2 + w_z/2, - b_gyro_y/2 - w_y/2, 0, 0, 0,  q0/2, -q3/2,  q2/2, 0];
            F(9, :) = [0, 0, 0, 0, 0, 0, b_gyro_y/2 + w_y/2, - b_gyro_z/2 - w_z/2,                 1/2,   b_gyro_x/2 + w_x/2, 0, 0, 0,  q3/2,  q0/2, -q1/2, 0];
            F(10, :) = [0, 0, 0, 0, 0, 0, b_gyro_z/2 + w_z/2,   b_gyro_y/2 + w_y/2, - b_gyro_x/2 - w_x/2,                 1/2, 0, 0, 0, -q2/2,  q1/2,  q0/2, 0];
        end

        function G = get_G_matrix(obj, x, u, w)
            G = zeros(17, 14);
            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);

            %dv_dot 
            G(4, :) = [          1 - 2*conj(q3)^2 - 2*conj(q2)^2, 2*conj(q1)*conj(q2) - 2*conj(q0)*conj(q3), 2*conj(q0)*conj(q2) + 2*conj(q1)*conj(q3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            G(5, :) = [2*conj(q0)*conj(q3) + 2*conj(q1)*conj(q2),           1 - 2*conj(q3)^2 - 2*conj(q1)^2, 2*conj(q2)*conj(q3) - 2*conj(q0)*conj(q1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            G(6, :) = [2*conj(q1)*conj(q3) - 2*conj(q0)*conj(q2), 2*conj(q0)*conj(q1) + 2*conj(q2)*conj(q3),           1 - 2*conj(q2)^2 - 2*conj(q1)^2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
           
            %d_q
            G(7, :) = [0, 0, 0, -q1/2, -q2/2, -q3/2, 0, 0, 0, 0, 0, 0, 0, 0];
            G(8, :) = [0, 0, 0,  q0/2, -q3/2,  q2/2, 0, 0, 0, 0, 0, 0, 0, 0];
            G(9, :) = [0, 0, 0,  q3/2,  q0/2, -q1/2, 0, 0, 0, 0, 0, 0, 0, 0];
            G(10, :) = [0, 0, 0, -q2/2,  q1/2,  q0/2, 0, 0, 0, 0, 0, 0, 0, 0];
            
            %dbeta^dot
            G(11:17, 8:14) = eye(7);
        end

        function z_hat  = get_measurement_estimate(obj, x)
            z_hat = x(3);
        end

        function H = get_H_matrix(obj)
            H = [...
                0, 0, 1,    0, 0, 0,    0, 0, 0, 0,    0, 0, 0    0, 0, 0,   0];
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            scale_var = 0.5*(1./(Fs.^2));
            vel_delta_bias_sigma = scale_var.* obj.accel_bias_noise;
            ang_vel_delta_bias_sigma = scale_var.* obj.gyro_bias_noise;
            pos_delta_bias_sigma = scale_var.* obj.baro_bias_noise;

            Qs = diag([obj.additive_noise.*ones(1,3), obj.additive_noise.*ones(1,3), obj.additive_noise.*ones(1,1), vel_delta_bias_sigma.*ones(1,3), ang_vel_delta_bias_sigma.*ones(1,3), pos_delta_bias_sigma.*ones(1,1)]);
            w = scale_var.*[obj.accel_noise.*ones(1,3), obj.gyro_noise.*ones(1,3), obj.baro_noise*ones(1,1), obj.accel_bias_noise.*ones(1,3), obj.gyro_bias_noise.*ones(1,3), obj.baro_bias_noise.*ones(1,1)];
        end

        function R = get_R_matrix(obj)
            R = obj.baro_measurement_uncertainty;
        end
    end
end

