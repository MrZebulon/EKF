classdef UKFEngine < BaseKFEngine
    properties (Constant)
        alpha = 1e-3;
        beta = 2;
        kappa = 1;
       
        sigma_x = 0;
    end
    methods
        function obj = predict_step(obj, u)
            N = obj.state_size();
        
            [Qs, w] = obj.model.generate_noise();
            Q = obj.model.get_Q_matrix(obj.x, u, w);
            [sigma, Wa, Wc] = obj.generate_sigma_points(N);

            % Apply f to each sigma point
            obj.sigma_x = arrayfun(@(i) obj.model.state_prediction_function(sigma(i, :), u), 1:2*N+1);

            % compute x(k+1,k) and P(k+1,k)
            obj.x = sum(Wa .* obj.sigma_x, 2);
            obj.P = cov_matrix(obj.sigma_x, obj.x, Wc, N) + Q + Qs;
        end

        function obj = update_step(obj, z, obs_func, R)
            N = obj.state_size();
            [sigma_z, Wa, Wc] = obj.generate_sigma_points(N);

            % Apply f to each sigma point
            sigma_z = arrayfun(@(i) obs_func(sigma_z(i, :), u), 1:2*N+1);
            
            % compute z_hat(k,k), S(k) and Cxz
            z_hat = sum(Wa .* sigma_z, 2);
            S = cov_matrix(sigma_z, z_hat, Wc, N) + R;
            Cxz = cov_matrixx(simga_z, z_hat, obj.sigma_x, obj.x, Wc, N);
            
            %Update
            K = Cxz*S;
            obj.x = obj.x + K*(z_hat-z);
            obj.P = obj.P + K*S*(K');
        end
        
        function obj = update_step_baro(obj, z)
            obj.update_step(z, obj.model.baro_obs_function, obj.model.baro_R_matrix(obj.x));
        end

        function obj = update_step_mag(obj, z)
             obj.update_step(z, obj.model.mag_obs_function, obj.model.mag_R_matrix(obj.x));
        end
        
        function Nx = state_size(obj)
            Nx = obj.model.Nx;
        end

        function [S, Wa, Wc] = generate_sigma_points(obj, N)
            S = zeros(2*N+1, 2*N+1);
            Wa = zeros(2*N+1, 1);
            Wc = zeros(2*N+1, 1);
            A = chol(obj.P, "lower");

            S(1, :) = (obj.x).';
            Wa(1) = (obj.alpha.^2 * obj.kappa - N)/(obj.alpha.^2*obj.kappa);
            Wc(1) = Wa(1) + 1 - obj.alpha.^2 + obj.beta;
            
            for j = 2:N
                S(j, :) = (obj.x + obj.alpha*sqrt(obj.kappa) * A(:, j)).';
                Wa(j) = 1/(2*obj.alpha.^2*obj.kappa);
                Wc(j) = 1/(2*obj.alpha.^2*obj.kappa);
                
                S(N+j, :) = (obj.x - obj.alpha*sqrt(obj.kappa) * A(:, j)).';
                Wa(N+j) = 1/(2*obj.alpha.^2*obj.kappa);
                Wc(N+j) = 1/(2*obj.alpha.^2*obj.kappa);
            end
        end

        function M = cov_matrix(sigma, mean, w, n)
            M = cov_matrixx(sigma, mean, sigma, mean, w, n);
        end

        function M = cov_matrixx(sigma1, mean1, sigma2, mean2, w, n)
            for i = 1:2*n+1
                M = M + w(i)*(sigma1(:, i) - mean1) * (sigma2(:, i) - mean2).';
            end
        end
    end
end
