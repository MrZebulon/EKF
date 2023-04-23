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
        
        % units = g (9.81 m/s)
        accel_bias = [-0.2668811881 0.3305544554 0.07695049505];
        accel_noise = 7.05E-04;
        accel_bias_noise =  6.89e-4;

        gyro_bias = [0 0 0];
        gyro_noise = 0;
        gyro_bias_noise = 0;

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

                0
                0
                0
                0
                
                obj.accel_bias(1)
                obj.accel_bias(2)
                obj.accel_bias(3)

                obj.gyro_bias(1)
                obj.gyro_bias(2)
                obj.gyro_bias(3)

                obj.baro_bias];

            P_init = diag(1e-9 * ones(1, 17));
        end

        function dx = get_delta_x(obj, x, u)
            dt = obj.dt;

            delta_a = Utils.body_to_inertial_frame(x(7:10), u(1:3) + x(11:16));
            delta_q = 0.5 * Utils.mult_quat(x(7:10), u(4:6) + x(14:16));
            
            dx = [
                x(4) * dt
                x(5) * dt
                x(6) * dt

                delta_a(1) * dt
                delta_a(2) * dt
                delta_a(3) * dt

                delta_q(1) * dt
                delta_q(2) * dt
                delta_q(3) * dt
                delta_q(4) * dt
                
                0
                0
                0

                0
                0
                0
                
                0];
        end

        function F = get_F_matrix(obj, x, u)
            
            F = zeros(17, 17);

            [q1, q2, q3, q4] = x(7:10);
            [a_x, a_y, a_z] = u(1:3);
            [w_x, w_y, w_z] = u(4:6);
            [b_acc_x, b_acc_y, b_acc_z] = x(11:13);
            [b_gyro_x, b_gyro_y, b_gyro_z] = x(14:16);

            % dp^dot/dv
            F(1:3, 4:6) = eye(3);
            
            %dv^dot/dq
            F(4, 7:10) = [2*q3*(a_z + b_acc_z) - 2*q4*(a_y + b_acc_y), 2*q3*(a_y + b_acc_y) + 2*q4*(a_z + b_acc_z), 2*q2*(a_y + b_acc_y) - 4*q3*(a_x + b_acc_x) + 2*q1*(a_z + b_acc_z), 2*q2*(a_z + b_acc_z) - 2*q1*(a_y + b_acc_y) - 4*q4*(a_x + b_acc_x)];
            F(5, 7:10) = [2*q4*(a_x + b_acc_x) - 2*q2*(a_z + b_acc_z), 2*q3*(a_x + b_acc_x) - 4*q2*(a_y + b_acc_y) - 2*q1*(a_z + b_acc_z), 2*q2*(a_x + b_acc_x) + 2*q4*(a_z + b_acc_z), 2*q1*(a_x + b_acc_x) - 4*q4*(a_y + b_acc_y) + 2*q3*(a_z + b_acc_z)];
            F(6, 7:10) = [2*q2*(a_y + b_acc_y) - 2*q3*(a_x + b_acc_x), 2*q4*(a_x + b_acc_x) + 2*q1*(a_y + b_acc_y) - 4*q2*(a_z + b_acc_z), 2*q4*(a_y + b_acc_y) - 2*q1*(a_x + b_acc_x) - 4*q3*(a_z + b_acc_z), 2*q2*(a_x + b_acc_x) + 2*q3*(a_y + b_acc_y)];
            
            %dv^dot/dbeta_acc
            F(4, 11:13) = [- 2*q3^2 - 2*q4^2 + 1, 2*q2*q3 - 2*q1*q4, 2*q1*q3 + 2*q2*q4];
            F(5, 11:13) = [2*q1*q4 + 2*q2*q3, - 2*q2^2 - 2*q4^2 + 1, 2*q3*q4 - 2*q1*q2];
            F(6, 11:13) = [2*q2*q4 - 2*q1*q3, 2*q1*q2 + 2*q3*q4, - 2*q2^2 - 2*q3^2 + 1];
            
            %dq^dot/dq
            F(7, 7:10) = [0, - b_gyro_x/2 - w_x/2, - b_gyro_y/2 - w_y/2, - b_gyro_z/2 - w_z/2];
            F(8, 7:10) = [b_gyro_x/2 + w_x/2, 0, b_gyro_z/2 + w_z/2, - b_gyro_y/2 - w_y/2];
            F(9, 7:10) = [b_gyro_y/2 + w_y/2, - b_gyro_z/2 - w_z/2, 0, b_gyro_x/2 + w_x/2];
            F(10, 7:10) = [b_gyro_z/2 + w_z/2, b_gyro_y/2 + w_y/2, - b_gyro_x/2 - w_x/2, 0];
            
            %dq^dot/dbeta_gyro
            F(7, 14:16) = [-q2/2, -q3/2, -q4/2];
            F(8, 14:16) = [q1/2, -q4/2,  q3/2];
            F(9, 14:16) = [q4/2,  q1/2, -q2/2];
            F(10, 14:16) = [-q3/2,  q2/2,  q1/2];
             
        end

        function G = get_G_matrix(obj, x, u, w)
            G = zeros(17, 14);
            
            %dv^dot/d_eta_acc
            G(4, 1:3) = [- 2*q3^2 - 2*q4^2 + 1, 2*q2*q3 - 2*q1*q4, 2*q1*q3 + 2*q2*q4];
            G(5, 1:3) = [2*q1*q4 + 2*q2*q3, - 2*q2^2 - 2*q4^2 + 1, 2*q3*q4 - 2*q1*q2];
            G(6, 1:3) = [2*q2*q4 - 2*q1*q3, 2*q1*q2 + 2*q3*q4, - 2*q2^2 - 2*q3^2 + 1];
            
            %dv^dot/d_eta_gyro
            G(7, 4:6) = [-q2/2, -q3/2, -q4/2];
            G(8, 4:6) = [q1/2, -q4/2,  q3/2];
            G(9, 4:6) = [q4/2,  q1/2, -q2/2];
            G(10, 4:6) = [-q3/2, q2/2, q1/2];
            
            %dbeta^dot/d_eta_beta
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

