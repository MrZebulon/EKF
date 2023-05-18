classdef TranslationRotationModelV2 < BaseModel
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
        
        % Time has been taken into account within the biases/noises/drifts
        % units = m/s
        accel_bias = [0.0029394893123127377, -0.0009818539510592773, 0.0028762637247315066];
        accel_noise = 4.7358607479753485e-09;
        accel_drift = 3.314312818032142e-10;

        %units = rad/s
        gyro_bias = [0.0038284331173227743, -0.001784190927555866, -0.002920534852288187];
        gyro_noise = 1.0102261028815603e-08;
        gyro_drift = 3.9979150848914756e-10;

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

                obj.baro_bias];

            P_init = diag(1e-9 * ones(1, 17));
        end

        function x_new = compute_x_new(obj, x, u)
            dt = obj.dt;
            delta_acc = Utils.body_to_inertial_frame(x(7:10)', (dt*u(1:3)' - x(11:13)));
            delta_q = quaternion([1, (dt * u(4:6) - x(14:16)')/ 2]);
            dx = [
                x(4) * dt
                x(5) * dt
                x(6) * dt

                delta_acc(1)
                delta_acc(2)
                delta_acc(3)

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
            x_new(7:10) = compact(normalize(quaternion(x(7:10)') * delta_q));
        end

        function F = get_F_matrix(obj, x, u)

            F = zeros(17, 17);

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            a_x = u(1) * obj.dt; a_y = u(2) * obj.dt; a_z = u(3) * obj.dt;
            w_x = u(4) * obj.dt; w_y = u(5) * obj.dt; w_z = u(6) * obj.dt;
            b_acc_x = x(11); b_acc_y = x(12); b_acc_z = x(13);
            b_gyro_x = x(14); b_gyro_y = x(15); b_gyro_z = x(16);

            % dp^dot
            F(1:3, 4:6) = eye(3) * obj.dt;
            
            %dv^dot
            F(4, :) = [0, 0, 0, 0, 0, 0, 2*q0*(a_x  - b_acc_x) - 2*q3*(a_y  - b_acc_y) + 2*q2*(a_z  - b_acc_z), 2*q1*(a_x  - b_acc_x) + 2*q2*(a_y  - b_acc_y) + 2*q3*(a_z  - b_acc_z), 2*q1*(a_y  - b_acc_y) - 2*q2*(a_x  - b_acc_x) + 2*q0*(a_z  - b_acc_z), 2*q1*(a_z  - b_acc_z) - 2*q0*(a_y  - b_acc_y) - 2*q3*(a_x  - b_acc_x), 0, 0, 0, - q0^2 - q1^2 + q2^2 + q3^2, 2*q0*q3 - 2*q1*q2, - 2*q0*q2 - 2*q1*q3, 0];
            F(5, :) = [0, 0, 0, 0, 0, 0, 2*q3*(a_x  - b_acc_x) + 2*q0*(a_y  - b_acc_y) - 2*q1*(a_z  - b_acc_z), 2*q2*(a_x  - b_acc_x) - 2*q1*(a_y  - b_acc_y) - 2*q0*(a_z  - b_acc_z), 2*q1*(a_x  - b_acc_x) + 2*q2*(a_y  - b_acc_y) + 2*q3*(a_z  - b_acc_z), 2*q0*(a_x  - b_acc_x) - 2*q3*(a_y  - b_acc_y) + 2*q2*(a_z  - b_acc_z), 0, 0, 0, - 2*q0*q3 - 2*q1*q2, - q0^2 + q1^2 - q2^2 + q3^2, 2*q0*q1 - 2*q2*q3,  0];
            F(6, :) = [0, 0, 0, 0, 0, 0, 2*q1*(a_y  - b_acc_y) - 2*q2*(a_x  - b_acc_x) + 2*q0*(a_z  - b_acc_z), 2*q3*(a_x  - b_acc_x) + 2*q0*(a_y  - b_acc_y) - 2*q1*(a_z  - b_acc_z), 2*q3*(a_y  - b_acc_y) - 2*q0*(a_x  - b_acc_x) - 2*q2*(a_z  - b_acc_z), 2*q1*(a_x  - b_acc_x) + 2*q2*(a_y  - b_acc_y) + 2*q3*(a_z  - b_acc_z), 0, 0, 0,   2*q0*q2 - 2*q1*q3, - 2*q0*q1 - 2*q2*q3, - q0^2 + q1^2 + q2^2 - q3^2, 0];
            
            %dq^dot
            F(7, :) = [0, 0, 0, 0, 0, 0, 0, b_gyro_x/2 - w_x/2, b_gyro_y/2 - w_y/2, b_gyro_z/2 - w_z/2, 0, 0, 0, q1/2, q2/2, q3/2, 0];
            F(8, :) = [0, 0, 0, 0, 0, 0, w_x/2 - b_gyro_x/2, 0, w_z/2 - b_gyro_z/2, b_gyro_y/2 - w_y/2, 0, 0, 0, -q0/2, q3/2, -q2/2, 0];
            F(9, :) = [0, 0, 0, 0, 0, 0, w_y/2 - b_gyro_y/2, b_gyro_z/2 - w_z/2, 0, w_x/2 - b_gyro_x/2, 0, 0, 0, -q3/2, -q0/2, q1/2, 0];
            F(10, :) = [0, 0, 0, 0, 0, 0, w_z/2 - b_gyro_z/2, w_y/2 - b_gyro_y/2, b_gyro_x/2 - w_x/2, 0, 0, 0, 0, q2/2, -q1/2, -q0/2, 0];
        end

        function Q = get_Q_matrix(obj, x, u, w)
            G = zeros(17, 7);
            G(4:6, 1:3) = Utils.rotmat(x(7:10));
            G(7:10, 4:7) = 0.5.*Utils.hamilton_product_as_matrix(x(7:10));

            Q = G*diag([w(1:3), 0, w(4:6)])*(G.');
        end

        function z_hat  = get_measurement_estimate(obj, x)
            % If using multiple sensors here, h just gets more columns
            z_hat = [x(3) + x(17)];
        end

        function H = get_H_matrix(obj)
            % If using multiple sensors here, H just gets more rows
            H = [...
                0, 0, 1,    0, 0, 0,    0, 0, 0, 0,    0, 0, 0    0, 0, 0,   0];
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            scale_var = 0.5*(1./(Fs.^2));
            vel_delta_bias_sigma = scale_var.* obj.accel_drift;
            ang_vel_delta_bias_sigma = scale_var.* obj.gyro_drift;
            pos_delta_bias_sigma = scale_var.* obj.baro_drift;

            Qs = diag([obj.additive_noise.*ones(1,10), vel_delta_bias_sigma.*ones(1,3), ang_vel_delta_bias_sigma.*ones(1,3), pos_delta_bias_sigma.*ones(1,1)]);
            w = scale_var.*[obj.accel_noise.*ones(1,3), obj.gyro_noise.*ones(1,3), obj.baro_noise.*ones(1,1)];
        end

        function R = get_R_matrix(obj)
            R = obj.baro_noise;
        end
    end
end
