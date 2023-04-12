classdef MEKF
    properties
        % x: state vector
        % P: covariance martix
        %
        % f: state transition model
        % gen_F: generator for df/dx (Jacobian Matrix w/ respect to state vector)
        % gen_G: generator for df/dw (Jacobian matrix w/ respect to noise vector)
        %
        % h: observation model
        % gen_H: generator for dh/dx (Jacobian Matrix w/ respect to state vector)
        %
        % noise_gen: generator for noise vector and process noise matrix
        % R: measurement uncertainty matrix
        x;
        P;
        model;
    end
    
    methods
        function obj = MEKF(x_init, P_init, model)
            obj.x = x_init;
            obj.P = P_init;
            obj.model = model;
        end

        function obj = predict_step(obj, u)
            % Prediction = x(k-1|k-1) -> x(k|k-1)
            % i.e. it computes the "a priori" estimate

            F = obj.model.get_F_matrix(obj.x, u);
            A = eye(size(F)) + F; % adding identity matrix to ensure that process noise isn't mistakenly removed
            x_new = obj.model.get_delta_x(obj.x, u) + obj.x; % x = delta_x + x_ref

            [Qs, w] = obj.model.generate_noise();
            G = obj.model.get_G_matrix(obj.x, u, w);
            P_new = A*obj.P*(A') + G*Qs*(G');

            obj.x = x_new;
            obj.P = P_new;
        end
        
        function obj = update_step(obj, z)
            % Update = x(k|k-1) -> x(k|k)
            % i.e. it computes the "a posteriori" estimate

            H = obj.model.get_H_matrix();
            nx = size(obj.x, 1);
            inov = z - obj.model.get_measurement_estimate(obj.x);
            S = H*obj.P*(H') + obj.model.get_R_matrix();
            K = obj.P*(H')*inv(S);

            obj.x = obj.x + K*inov;
            obj.P = (eye(nx)-K*H)*obj.P;
        end

        function n = state_size(obj)
            n = size(obj.x, 1); % state is a n-by-1 array
        end
    end
end
