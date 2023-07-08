classdef TranslationRotationModelV2 < BaseModel
    %
    % State vector definition
    %   1-3 : position
    %   4-6 : velocity
    %   7-10 : orientation
    %   11-13 : accel bias
    %   14-16 : gyro bias
    %   17-19: magneto_bias
    %   20-22 : earth permanent field in inertial frame
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
        Nz = 4;
    end
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

        % units = uT
        magneto_bias = [100, 100, 100];
        magneto_noise = 0.09;
        magneto_drift =  1e-10;
        magneto_earth_field = [19.5281 -5.0741 48.0067];
        magneto_earth_field_drift =  1e-15;

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

            P_init = diag(1e-9 * ones(1, 22));
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
            delta_q = [1.0; (dt.*u(4:6)'-x(14:16))/2];
            x_new(7:10) = Utils.mult_quat(x(7:10),delta_q);
        end

        function F = get_F_matrix(obj, x, u)
            F = zeros(obj.Nx, obj.Nx);

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);

            dvx_b = x(11); dvy_b = x(12); dvz_b = x(13);
            dvx = u(1) * obj.dt; dvy = u(2) * obj.dt; dvz = u(3) * obj.dt;


            dax = u(4) * obj.dt; day = u(5) * obj.dt; daz = u(6) * obj.dt;
            dax_b = x(14); day_b = x(15); daz_b = x(16);

            F(1:3, 4:6) = obj.dt.*eye(3);

            F(4, 7:10) = [2*q0*(dvx - dvx_b) - 2*q3*(dvy - dvy_b) + 2*q2*(dvz - dvz_b), 2*q1*(dvx - dvx_b) + 2*q2*(dvy - dvy_b) + 2*q3*(dvz - dvz_b), 2*q1*(dvy - dvy_b) - 2*q2*(dvx - dvx_b) + 2*q0*(dvz - dvz_b), 2*q1*(dvz - dvz_b) - 2*q0*(dvy - dvy_b) - 2*q3*(dvx - dvx_b)];
            F(4, 14:16) = [- q0^2 - q1^2 + q2^2 + q3^2, 2*q0*q3 - 2*q1*q2, - 2*q0*q2 - 2*q1*q3];
            F(5, 7:10) = [2*q3*(dvx - dvx_b) + 2*q0*(dvy - dvy_b) - 2*q1*(dvz - dvz_b), 2*q2*(dvx - dvx_b) - 2*q1*(dvy - dvy_b) - 2*q0*(dvz - dvz_b), 2*q1*(dvx - dvx_b) + 2*q2*(dvy - dvy_b) + 2*q3*(dvz - dvz_b), 2*q0*(dvx - dvx_b) - 2*q3*(dvy - dvy_b) + 2*q2*(dvz - dvz_b)];
            F(5, 14:16) = [- 2*q0*q3 - 2*q1*q2, - q0^2 + q1^2 - q2^2 + q3^2, 2*q0*q1 - 2*q2*q3];
            F(6, 7:10) = [2*q1*(dvy - dvy_b) - 2*q2*(dvx - dvx_b) + 2*q0*(dvz - dvz_b), 2*q3*(dvx - dvx_b) + 2*q0*(dvy - dvy_b) - 2*q1*(dvz - dvz_b), 2*q3*(dvy - dvy_b) - 2*q0*(dvx - dvx_b) - 2*q2*(dvz - dvz_b), 2*q1*(dvx - dvx_b) + 2*q2*(dvy - dvy_b) + 2*q3*(dvz - dvz_b)];
            F(6, 14:16) = [2*q0*q2 - 2*q1*q3,  - 2*q0*q1 - 2*q2*q3, - q0^2 + q1^2 + q2^2 - q3^2];

            F(7, 7:10) = [0,  dax_b/2 - dax/2,  day_b/2 - day/2,  daz_b/2 - daz/2];
            F(7, 14:16) = [q1/2,  q2/2,  q3/2];
            F(8, 7:10) = [dax/2 - dax_b/2, 0,  daz/2 - daz_b/2,  day_b/2 - day/2];
            F(8, 14:16) = [-q0/2,  q3/2, -q2/2];
            F(9, 7:10) = [day/2 - day_b/2,  daz_b/2 - daz/2, 0,  dax/2 - dax_b/2];
            F(9, 14:16) = [-q3/2, -q0/2,  q1/2];
            F(10, 7:10) = [daz/2 - daz_b/2,  day/2 - day_b/2,  dax_b/2 - dax/2, 0];
            F(10, 7:10) = [q2/2, -q1/2, -q0/2];
        end

        function Q = get_Q_matrix(obj, x, u, w)

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            G = zeros(obj.Nx, obj.Nw);

            G(4:6, 1:3) = quat2rotm([q0, q1, q2, q3]);
            G(7:10, 4:7) = 0.5 .*Utils.hamilton_product_as_matrix([q0, q1, q2, q3]);
            Q = G*diag([w(1:3), w(4:7)])*(G.');

        end

        function z_hat  = get_measurement_estimate(obj, x)

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            magNavX = obj.x(17); magNavY = obj.x(18); magNavZ = obj.x(19);
            magBiasX = obj.x(20); magBiasY = obj.x(21); magBiasZ = obj.x(22);

            mx = magBiasX + magNavX*(q0^2 + q1^2 - q2^2 - q3^2) - magNavZ*(2*q0*q2 - 2*q1*q3) + magNavY*(2*q0*q3 + 2*q1*q2);
            my = magBiasY + magNavY*(q0^2 - q1^2 + q2^2 - q3^2) + magNavZ*(2*q0*q1 + 2*q2*q3) - magNavX*(2*q0*q3 - 2*q1*q2);
            mz = magBiasZ + magNavZ*(q0^2 - q1^2 - q2^2 + q3^2) - magNavY*(2*q0*q1 - 2*q2*q3) + magNavX*(2*q0*q2 + 2*q1*q3);

            z_hat = [x(3) + obj.baro_bias, mx, my, mz];
        end

        function H = get_H_matrix(obj)
            % If using multiple sensors here, H just gets more rows
            H = Zeros(obj.Nz, obj.Nx);

            q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
            magNavX = obj.x(17); magNavY = obj.x(18); magNavZ = obj.x(19);

            H(1, 3) = 1.0;

            H(2, 7:10) = [2*magNavY*q3 - 2*magNavZ*q2 + 2*magNavX*q0, 2*magNavZ*q3 + 2*magNavY*q2 + 2*magNavX*q1, 2*magNavY*q1 - 2*magNavZ*q0 - 2*magNavX*q2, 2*magNavZ*q1 + 2*magNavY*q0 - 2*magNavX*q3];
            H(2, 17:19) = [q0^2 + q1^2 - q2^2 - q3^2, 2*q0*q3 + 2*q1*q2, 2*q1*q3 - 2*q0*q2];
            H(3, 7:10) = [2*magNavZ*q1 + 2*magNavY*q0 - 2*magNavX*q3, 2*magNavZ*q0 - 2*magNavY*q1 + 2*magNavX*q2, 2*magNavZ*q3 + 2*magNavY*q2 + 2*magNavX*q1, 2*magNavZ*q2 - 2*magNavY*q3 - 2*magNavX*q0];
            H(3, 17:19) = [2*q1*q2 - 2*q0*q3, q0^2 - q1^2 + q2^2 - q3^2, 2*q0*q1 + 2*q2*q3];
            H(4, 7:10) = [2*magNavZ*q0 - 2*magNavY*q1 + 2*magNavX*q2, 2*magNavX*q3 - 2*magNavY*q0 - 2*magNavZ*q1, 2*magNavY*q3 - 2*magNavZ*q2 + 2*magNavX*q0, 2*magNavZ*q3 + 2*magNavY*q2 + 2*magNavX*q1];
            H(4, 17:19) = [2*q0*q2 + 2*q1*q3, 2*q2*q3 - 2*q0*q1, q0^2 - q1^2 - q2^2 + q3^2];
            H(2:4, 20:22) = eye(3);
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            scale_var = 0.5*(1./(Fs.^2));
            accel_drift_sigma = scale_var.* obj.accel_drift;
            gyro_drift_sigma = scale_var.* obj.gyro_drift;

            Qs = diag([obj.additive_noise.*ones(1,10), accel_drift_sigma.*ones(1,3), gyro_drift_sigma.*ones(1,3), obj.magneto_earth_field_drift.*ones(1,3), obj.magneto_drift.*ones(1,3)]);

            w = scale_var.*[obj.accel_noise.*ones(1,3), 0, obj.gyro_noise.*ones(1,3)];
        end

        function R = get_R_matrix(obj)
            R = diag([obj.baro_noise, obj.magneto_noise.*ones(1,3)]);
        end

    end
end
