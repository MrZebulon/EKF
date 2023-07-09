classdef BaseKFEngine < handle
   properties
        % x: state vector
        % P: covariance martix
        x;
        P;
        model;
    end
    methods
        function obj = BaseKFEngine(model)
            [x_init, P_init] = model.get_init_state();
            obj.x = x_init;
            obj.P = P_init;
            obj.model = model;
        end

        function obj = predict_step(obj)
            % Prediction = x(k-1|k-1) -> x(k|k-1)
            % i.e. it computes the "a priori" estimate
        end

        function obj = update_step(obj)
            % Update = x(k|k-1) -> x(k|k)
            % i.e. it computes the "a posteriori" estimate
        end
    end
end