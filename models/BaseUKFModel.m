classdef BaseUKFModel < handle
    properties
        dt;
    end
    methods
        function obj = BaseUKFModel(dt)
            obj.dt = dt;
        end
        
        function [x_init, P_init] = get_init_state(obj)
            x_init = 0;
            P_init = 0;
        end
        
        function x_new = state_prediction_function(obj, x, u)
            x_new = x;
        end

        function Q = get_Q_matrix(obj, x, u, w)
            Q = 0;
        end

        function [Qs, w] = generate_noise(obj)
            Qs = 0;
            w = 0;
        end

        function h = baro_obs_function(obj, x)
            h = 0;
        end

        function R = baro_R_matrix(obj, x)
            R = 0;
        end

        function h = mag_obs_function(obj, x)
            h = 0;
        end

        function R = mag_R_matrix(obj, x)
            R = 0;
        end

    end

end