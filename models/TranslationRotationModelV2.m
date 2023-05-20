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
        additive_noise = 1e-6;

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
                obj.gyro_bias(3)];

            P_init = diag(1e-9 * ones(1, 16));
        end

        function x_new = compute_x_new(obj, x, u)
            dt = obj.dt;
            delta_acc = Utils.rotmat(x(7:10)') * (dt*u(1:3)' - x(11:13));
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
                0];

            x_new = x + dx;

            % delta_q = quaternion((dt * u(4:6) - x(14:16)'), "rotvec");
            %x_new(7:10) = compact(normalize(quaternion(x(7:10)') * delta_q));
            delta_q = [1; (dt .* u(4:6)'-x(14:16))/2];
            x_new(7:10) = Utils.mult_quat(x(7:10),delta_q);
        end

        function F = get_F_matrix(obj, x, u)
            F = zeros(16, 16);

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);

            dvx_b = x(11); dvy_b = x(12); dvz_b = x(13);
            dvx = u(1) * obj.dt; dvy = u(2) * obj.dt; dvz = u(3) * obj.dt;


            dax = u(4) * obj.dt; day = u(5) * obj.dt; daz = u(6) * obj.dt;
            dax_b = x(14); day_b = x(15); daz_b = x(16);

            F(1:3, 4:6) = obj.dt.*eye(3);

            F(4, :) = [0, 0, 0, 0, 0, 0, 2*q0*(dvx - dvx_b) - 2*q3*(dvy - dvy_b) + 2*q2*(dvz - dvz_b), 2*q1*(dvx - dvx_b) + 2*q2*(dvy - dvy_b) + 2*q3*(dvz - dvz_b), 2*q1*(dvy - dvy_b) - 2*q2*(dvx - dvx_b) + 2*q0*(dvz - dvz_b), 2*q1*(dvz - dvz_b) - 2*q0*(dvy - dvy_b) - 2*q3*(dvx - dvx_b), 0, 0, 0, - q0^2 - q1^2 + q2^2 + q3^2, 2*q0*q3 - 2*q1*q2, - 2*q0*q2 - 2*q1*q3];
            F(5, :) = [0, 0, 0, 0, 0, 0, 2*q3*(dvx - dvx_b) + 2*q0*(dvy - dvy_b) - 2*q1*(dvz - dvz_b), 2*q2*(dvx - dvx_b) - 2*q1*(dvy - dvy_b) - 2*q0*(dvz - dvz_b), 2*q1*(dvx - dvx_b) + 2*q2*(dvy - dvy_b) + 2*q3*(dvz - dvz_b), 2*q0*(dvx - dvx_b) - 2*q3*(dvy - dvy_b) + 2*q2*(dvz - dvz_b), 0, 0, 0,  - 2*q0*q3 - 2*q1*q2, - q0^2 + q1^2 - q2^2 + q3^2, 2*q0*q1 - 2*q2*q3];
            F(6, :) = [0, 0, 0, 0, 0, 0, 2*q1*(dvy - dvy_b) - 2*q2*(dvx - dvx_b) + 2*q0*(dvz - dvz_b), 2*q3*(dvx - dvx_b) + 2*q0*(dvy - dvy_b) - 2*q1*(dvz - dvz_b), 2*q3*(dvy - dvy_b) - 2*q0*(dvx - dvx_b) - 2*q2*(dvz - dvz_b), 2*q1*(dvx - dvx_b) + 2*q2*(dvy - dvy_b) + 2*q3*(dvz - dvz_b), 0, 0, 0, 2*q0*q2 - 2*q1*q3,  - 2*q0*q1 - 2*q2*q3, - q0^2 + q1^2 + q2^2 - q3^2];

            F(7, :) = [0, 0, 0, 0, 0, 0, 0,  dax_b/2 - dax/2,  day_b/2 - day/2,  daz_b/2 - daz/2, 0, 0, 0,  q1/2,  q2/2,  q3/2];
            F(8, :) = [0, 0, 0, 0, 0, 0, dax/2 - dax_b/2, 0,  daz/2 - daz_b/2,  day_b/2 - day/2, 0, 0, 0, -q0/2,  q3/2, -q2/2];
            F(9, :) = [0, 0, 0, 0, 0, 0, day/2 - day_b/2,  daz_b/2 - daz/2, 0,  dax/2 - dax_b/2, 0, 0, 0, -q3/2, -q0/2,  q1/2];
            F(10, :) = [0, 0, 0, 0, 0, 0, daz/2 - daz_b/2,  day/2 - day_b/2,  dax_b/2 - dax/2, 0, 0, 0, 0,  q2/2, -q1/2, -q0/2];
        end

        function Q = get_Q_matrix(obj, x, u, w)

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            G = zeros(16, 7);

            G(4:6, 1:3) = quat2rotm([q0, q1, q2, q3]);
            G(7:10, 4:7) = 0.5 .*Utils.hamilton_product_as_matrix([q0, q1, q2, q3]);
            Q = G*diag([w(1:3), 0, w(4:6)])*(G.');

        end
        
        %{
        function F = get_F_matrixAJD(obj, x, u)

            F = zeros(16, 16);

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            a_x = u(1) * obj.dt; a_y = u(2) * obj.dt; a_z = u(3) * obj.dt;
            w_x = u(4) * obj.dt; w_y = u(5) * obj.dt; w_z = u(6) * obj.dt;
            b_acc_x = x(11); b_acc_y = x(12); b_acc_z = x(13);
            b_gyro_x = x(14); b_gyro_y = x(15); b_gyro_z = x(16);

            % dp^dot
            F(1:3, 4:6) = eye(3) * obj.dt;

            %dv^dot
            F(4, :) = [0, 0, 0, 0, 0, 0, 2*q0*(a_x  - b_acc_x) - 2*q3*(a_y  - b_acc_y) + 2*q2*(a_z  - b_acc_z), 2*q1*(a_x  - b_acc_x) + 2*q2*(a_y  - b_acc_y) + 2*q3*(a_z  - b_acc_z), 2*q1*(a_y  - b_acc_y) - 2*q2*(a_x  - b_acc_x) + 2*q0*(a_z  - b_acc_z), 2*q1*(a_z  - b_acc_z) - 2*q0*(a_y  - b_acc_y) - 2*q3*(a_x  - b_acc_x), 0, 0, 0, - q0^2 - q1^2 + q2^2 + q3^2, 2*q0*q3 - 2*q1*q2, - 2*q0*q2 - 2*q1*q3];
            F(5, :) = [0, 0, 0, 0, 0, 0, 2*q3*(a_x  - b_acc_x) + 2*q0*(a_y  - b_acc_y) - 2*q1*(a_z  - b_acc_z), 2*q2*(a_x  - b_acc_x) - 2*q1*(a_y  - b_acc_y) - 2*q0*(a_z  - b_acc_z), 2*q1*(a_x  - b_acc_x) + 2*q2*(a_y  - b_acc_y) + 2*q3*(a_z  - b_acc_z), 2*q0*(a_x  - b_acc_x) - 2*q3*(a_y  - b_acc_y) + 2*q2*(a_z  - b_acc_z), 0, 0, 0, - 2*q0*q3 - 2*q1*q2, - q0^2 + q1^2 - q2^2 + q3^2, 2*q0*q1 - 2*q2*q3];
            F(6, :) = [0, 0, 0, 0, 0, 0, 2*q1*(a_y  - b_acc_y) - 2*q2*(a_x  - b_acc_x) + 2*q0*(a_z  - b_acc_z), 2*q3*(a_x  - b_acc_x) + 2*q0*(a_y  - b_acc_y) - 2*q1*(a_z  - b_acc_z), 2*q3*(a_y  - b_acc_y) - 2*q0*(a_x  - b_acc_x) - 2*q2*(a_z  - b_acc_z), 2*q1*(a_x  - b_acc_x) + 2*q2*(a_y  - b_acc_y) + 2*q3*(a_z  - b_acc_z), 0, 0, 0,   2*q0*q2 - 2*q1*q3, - 2*q0*q1 - 2*q2*q3, - q0^2 + q1^2 + q2^2 - q3^2];

            %dq^dot
            F(7, :) = [0, 0, 0, 0, 0, 0, 0, b_gyro_x/2 - w_x/2, b_gyro_y/2 - w_y/2, b_gyro_z/2 - w_z/2, 0, 0, 0, q1/2, q2/2, q3/2];
            F(8, :) = [0, 0, 0, 0, 0, 0, w_x/2 - b_gyro_x/2, 0, w_z/2 - b_gyro_z/2, b_gyro_y/2 - w_y/2, 0, 0, 0, -q0/2, q3/2, -q2/2];
            F(9, :) = [0, 0, 0, 0, 0, 0, w_y/2 - b_gyro_y/2, b_gyro_z/2 - w_z/2, 0, w_x/2 - b_gyro_x/2, 0, 0, 0, -q3/2, -q0/2, q1/2];
            F(10, :) = [0, 0, 0, 0, 0, 0, w_z/2 - b_gyro_z/2, w_y/2 - b_gyro_y/2, b_gyro_x/2 - w_x/2, 0, 0, 0, 0, q2/2, -q1/2, -q0/2];
        end
       
        function Q = get_Q_matrix(obj, x, u, w)
            Q = zeros(16, 16);

            dvxCov = w(1);
            dvyCov = w(2);
            dvzCov = w(3);

            daxCov = w(4);
            dayCov = w(5);
            dazCov = w(6);


            q0 = x(7);
            q1 = x(8);
            q2 = x(9);
            q3 = x(10);

            Q(7:10, 7:10) = [...
                (daxCov*q1^2)/4 + (dayCov*q2^2)/4 + (dazCov*q3^2)/4, (dayCov*q2*q3)/4 - (daxCov*q0*q1)/4 - (dazCov*q2*q3)/4, (dazCov*q1*q3)/4 - (dayCov*q0*q2)/4 - (daxCov*q1*q3)/4, (daxCov*q1*q2)/4 - (dayCov*q1*q2)/4 - (dazCov*q0*q3)/4
                (dayCov*q2*q3)/4 - (daxCov*q0*q1)/4 - (dazCov*q2*q3)/4,    (daxCov*q0^2)/4 + (dazCov*q2^2)/4 + (dayCov*q3^2)/4, (daxCov*q0*q3)/4 - (dayCov*q0*q3)/4 - (dazCov*q1*q2)/4, (dazCov*q0*q2)/4 - (dayCov*q1*q3)/4 - (daxCov*q0*q2)/4
                (dazCov*q1*q3)/4 - (dayCov*q0*q2)/4 - (daxCov*q1*q3)/4, (daxCov*q0*q3)/4 - (dayCov*q0*q3)/4 - (dazCov*q1*q2)/4,    (dayCov*q0^2)/4 + (dazCov*q1^2)/4 + (daxCov*q3^2)/4, (dayCov*q0*q1)/4 - (daxCov*q2*q3)/4 - (dazCov*q0*q1)/4
                (daxCov*q1*q2)/4 - (dayCov*q1*q2)/4 - (dazCov*q0*q3)/4, (dazCov*q0*q2)/4 - (dayCov*q1*q3)/4 - (daxCov*q0*q2)/4, (dayCov*q0*q1)/4 - (daxCov*q2*q3)/4 - (dazCov*q0*q1)/4,    (dazCov*q0^2)/4 + (dayCov*q1^2)/4 + (daxCov*q2^2)/4];

            Q(4:6, 4:6) = [...
                dvyCov*(2*q0*q3 - 2*q1*q2)^2 + dvzCov*(2*q0*q2 + 2*q1*q3)^2 + dvxCov*(q0^2 + q1^2 - q2^2 - q3^2)^2, dvxCov*(2*q0*q3 + 2*q1*q2)*(q0^2 + q1^2 - q2^2 - q3^2) - dvyCov*(2*q0*q3 - 2*q1*q2)*(q0^2 - q1^2 + q2^2 - q3^2) - dvzCov*(2*q0*q1 - 2*q2*q3)*(2*q0*q2 + 2*q1*q3), dvzCov*(2*q0*q2 + 2*q1*q3)*(q0^2 - q1^2 - q2^2 + q3^2) - dvxCov*(2*q0*q2 - 2*q1*q3)*(q0^2 + q1^2 - q2^2 - q3^2) - dvyCov*(2*q0*q1 + 2*q2*q3)*(2*q0*q3 - 2*q1*q2)
                dvxCov*(2*q0*q3 + 2*q1*q2)*(q0^2 + q1^2 - q2^2 - q3^2) - dvyCov*(2*q0*q3 - 2*q1*q2)*(q0^2 - q1^2 + q2^2 - q3^2) - dvzCov*(2*q0*q1 - 2*q2*q3)*(2*q0*q2 + 2*q1*q3), dvxCov*(2*q0*q3 + 2*q1*q2)^2 + dvzCov*(2*q0*q1 - 2*q2*q3)^2 + dvyCov*(q0^2 - q1^2 + q2^2 - q3^2)^2, dvyCov*(2*q0*q1 + 2*q2*q3)*(q0^2 - q1^2 + q2^2 - q3^2) - dvzCov*(2*q0*q1 - 2*q2*q3)*(q0^2 - q1^2 - q2^2 + q3^2) - dvxCov*(2*q0*q2 - 2*q1*q3)*(2*q0*q3 + 2*q1*q2)
                dvzCov*(2*q0*q2 + 2*q1*q3)*(q0^2 - q1^2 - q2^2 + q3^2) - dvxCov*(2*q0*q2 - 2*q1*q3)*(q0^2 + q1^2 - q2^2 - q3^2) - dvyCov*(2*q0*q1 + 2*q2*q3)*(2*q0*q3 - 2*q1*q2), dvyCov*(2*q0*q1 + 2*q2*q3)*(q0^2 - q1^2 + q2^2 - q3^2) - dvzCov*(2*q0*q1 - 2*q2*q3)*(q0^2 - q1^2 - q2^2 + q3^2) - dvxCov*(2*q0*q2 - 2*q1*q3)*(2*q0*q3 + 2*q1*q2), dvxCov*(2*q0*q2 - 2*q1*q3)^2 + dvyCov*(2*q0*q1 + 2*q2*q3)^2 + dvzCov*(q0^2 - q1^2 - q2^2 + q3^2)^2];
        end
        %}

        function z_hat  = get_measurement_estimate(obj, x)
            % If using multiple sensors here, h just gets more columns
            z_hat = [x(3) + obj.baro_bias];
        end

        function H = get_H_matrix(obj)
            % If using multiple sensors here, H just gets more rows
            H = [...
                0, 0, 1,    0, 0, 0,    0, 0, 0, 0,    0, 0, 0    0, 0, 0];
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            scale_var = 0.5*(1./(Fs.^2));
            vel_delta_bias_sigma = scale_var.* obj.accel_drift;
            ang_vel_delta_bias_sigma = scale_var.* obj.gyro_drift;

            Qs = diag([obj.additive_noise.*ones(1,10), vel_delta_bias_sigma.*ones(1,3), ang_vel_delta_bias_sigma.*ones(1,3)]);

            w = scale_var.*[obj.accel_noise.*ones(1,3), obj.gyro_noise.*ones(1,3)];
        end

        function R = get_R_matrix(obj)
            R = 0.1; %obj.baro_noise;
        end

        function x = callback(obj, engine)
            %anti-drift
            x = engine.x;
            x(14:16) = obj.gyro_bias;
        end
    end
end
