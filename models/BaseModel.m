classdef Model < handle
    properties
        dt;
    end
    methods
        function obj = Model(dt)
            obj.dt = dt;
        end
        function x_new = get_next_state(obj, x, u)
            % Discrete form of EKF : gives x(k|k-1)
            % It does not represent the f function in itself, it should
            % rather be seen as the following
            % x(k|k-1) = x(k-1|k-1) + x_new(k)
            %
            x_new = 0;
        end

        function F = get_F_matrix(obj, x, u)
            F = 0;
        end

        function G = get_G_matrix(obj, x, u, w)
            G = 0;
        end

        function z_hat  = get_measurement_estimate(obj, x)
            z_hat = 0;
        end

        function H = get_H_matrix(obj)
            H = 0;
        end

        function [Qs, w] = generate_noise(obj)
            Qs = 0;
            w = 0;
        end

        function R = get_R_matrix(obj)
            R = 0;
        end
    end

end