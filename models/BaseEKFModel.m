classdef BaseEKFModel < handle
    properties
        dt;
    end
    methods
        function obj = BaseEKFModel(dt)
            obj.dt = dt;
        end
        
        function [x_init, P_init] = get_init_state(obj)
            x_init = 0;
            P_init = 0;
        end
        function x_new = compute_x_new(obj, x, u)
            % x_new = x + delta_x
            % delta_x is defined as the integral over a given timestep (dt) of
            % the state update function f(x, u, w).
            x_new = x;
        end

        function F = get_F_matrix(obj, x, u)
            % The F matrix is the Jacobian matrix of the state update
            % function f(x, u, w) with respect to the state vector x.
            F = 0;
        end

        function Q = get_Q_matrix(obj, x, u, w)
            % The Q matrix is the process Covariance Matrix.
            % Q(i,j) = Cov(x_i, x_j)
            Q = 0;
        end

        function [Qs, w] = generate_noise(obj)
            % Qs represent the process noise matrix.
            % w represents the process noise vector.
            % Both shall be generated at each step so as to allow noise to
            % change over time (i.e. w(k) would depend on w(k-1)).
            Qs = 0;
            w = 0;
        end

        function [h, H, R] = baro_obs_model(obj, x)
            % z_hat is defined as the estimate of some or all elements of
            % the state vector. This estimate is formed from the current state vector.
            % H is the Jacobian matrix of the obeservation
            % function h(x) with respect to the state vector x.
            % R is the process noise matrix 
            h = 0;
            H = 0;
            R = 0;
        end

        function [h, H, R] = mag_obs_model(obj, x)
            % z_hat is defined as the estimate of some or all elements of
            % the state vector. This estimate is formed from the current state vector.
            % H is the Jacobian matrix of the obeservation
            % function h(x) with respect to the state vector x.
            % R is the process noise matrix 
            h = 0;
            H = 0;
            R = 0;
        end
    end

end