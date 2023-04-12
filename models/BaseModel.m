classdef BaseModel < handle
    properties
        dt;
    end
    methods
        function obj = BaseModel(dt)
            obj.dt = dt;
        end
        function dx = get_delta_x(obj, x, u)
            % delta_x is defined as the integral over a given timestep (dt) of
            % the state update function f(x, u, w).
            dx = 0;
        end

        function F = get_F_matrix(obj, x, u)
            % The F matrix is the Jacobian matrix of the state update
            % function f(x, u, w) with respect to the state vector x.
            F = 0;
        end

        function G = get_G_matrix(obj, x, u, w)
            % The G matrix is the Jacobian matrix of the state update
            % function f(x, u, w) with respect to the noise vector w.
            G = 0;
        end

        function z_hat  = get_measurement_estimate(obj, x)
            % z_hat is defined as the estimate of some or all elements of
            % the state vector. This estimate is formed from the current
            % state vector.
            z_hat = 0;
        end

        function H = get_H_matrix(obj)
            % The G matrix is the Jacobian matrix of the state update
            % function h(x) with respect to the state vector x.
            H = 0;
        end

        function [Qs, w] = generate_noise(obj)
            % Qs represent the process noise matrix.
            % w represents the process noise vector.
            % Both shall be generated at each step so as to allow noise to
            % change over time (i.e. w(k) would depend on w(k-1)).
            Qs = 0;
            w = 0;
        end

        function R = get_R_matrix(obj)
            % R reprensents the measurement uncertainty matrix.
            % Generated at each step so as to allow it to
            % change over time (i.e. R(k) would depend on R(k-1)).
            R = 0;
        end
    end

end