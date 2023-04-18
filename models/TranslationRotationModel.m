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

    properties
        additive_noise = 1e-8;
        
        accel_bias = [0 0 0];
        accel_noise = 7.05E-04;
        accel_bias_noise =  6.89e-4;

        gyro_bias = [0 0 0];
        gyro_noise = 7.05E-04;
        gyro_bias_noise =  6.89e-4;

        baro_bias = 0;
        baro_noise = 1.52e-5; % units = hPa
        baro_bias_noise = 2.98e-7; % units = hPa
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
            F(1:3, 4:6) = eye(3);

            F(4, 7:10) = [-2*x(10)*(u(2) + x(12)) - 2*x(9)*(u(3) + x(13)), 2*x(9)*(u(2) + x(12)) + 2*x(10)*(u(3) + x(13)), 2*x(8)*(u(2) + x(12)) - 4*x(9)*(u(1) + x(11)) - 2*x(7)*(u(3) + x(13)), 2*x(8)*(u(3) + x(13)) - 2*x(7)*(u(2) + x(12)) - 4*x(10)*(u(1) + x(11))]';
            F(5, 7:10) = [ 2*x(10)*(u(1) + x(11)) - 2*x(8)*(u(3) + x(13)), 2*x(9)*(u(1) + x(11)) - 4*x(8)*(u(2) + x(12)) - 2*x(7)*(u(3) + x(13)), 2*x(8)*(u(1) + x(11)) + 2*x(10)*(u(3) + x(13)), 2*x(7)*(u(1) + x(11)) - 4*x(10)*(u(2) + x(12)) + 2*x(9)*(u(3) + x(13))]';
            F(6, 7:10) = [ 2*x(8)*(u(3) + x(13)) - 2*x(9)*(u(3) + x(13)), 2*x(10)*(u(3) + x(13)) + 2*x(7)*(u(3) + x(13)) - 4*x(8)*(u(3) + x(13)), 2*x(10)*(u(3) + x(13)) - 2*x(7)*(u(3) + x(13)) - 4*x(9)*(u(3) + x(13)), 2*x(8)*(u(3) + x(13)) + 2*x(9)*(u(3) + x(13))]';
         
            F(4, 11:13) = [ -2*x(9)^2 - 2*x(10)^2 + 1, 2*x(8)*x(9) - 2*x(7)*x(10), 2*x(8)*x(10) - 2*x(7)*x(9)]';
            F(5, 11:13) = [ 2*x(7)*x(10) + 2*x(8)*x(9), - 2*x(8)^2 - 2*x(10)^2 + 1, 2*x(9)*x(10) - 2*x(7)*x(8)]';
            F(6, 11:13) = [ 2*x(8)*x(10) - 2*x(7)*x(9), 2*x(7)*x(8) + 2*x(9)*x(10), - 2*x(8)^2 - 2*x(9)^2 + 1]';
           
            F(7, 7:10) = [ 0, (u(4) + x(14)), -(u(5) + x(15)), -(u(6) + x(16))]';
            F(8, 7:10) = [ (u(4) + x(14)), 0, (u(6) + x(16)), -(u(5) + x(15))]';
            F(9, 7:10) = [ (u(5) + x(15)), -(u(6) + x(16)), 0, (u(4) + x(14))]';
            F(10, 7:10) = [ (u(6) + x(16)), (u(5) + x(15)), -(u(4) + x(14)), 0]';
             
            F(14:16, 7) = [ x(8), -x(9), -x(10)]';
            F(14:16, 8) = [ x(7), -x(10),  x(9)]';
            F(14:16, 9) = [ x(10),  x(7), -x(8)]';
            F(14:16, 10) = [-x(9),  x(8),  x(7)]';

        end

        function G = get_G_matrix(obj, x, u, w)
            G = zeros(17, 14);
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

