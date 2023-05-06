classdef Utils
    methods (Static)
        function vect_i = body_to_inertial_frame(q, vect_b)
            rot_matrix = quat2rotm(q/norm(q));
            vect_i = rot_matrix * vect_b;
        end

        function qn = mult_quat(q1,q2)
            qnew_0 = q1(1)*q2(1)    -q1(2)*q2(2)    -q1(3)*q2(3)    -q1(4)*q2(4);
            qnew_1 = q1(1)*q2(2)    +q2(1)*q1(2)    +q1(3)*q2(4)    -q2(3)*q1(4);
            qnew_2 = q1(1)*q2(3)    +q2(1)*q1(3)    -q1(2)*q2(4)    +q2(2)*q1(4);
            qnew_3 = q1(1)*q2(4)    +q2(1)*q1(4)    +q1(2)*q2(3)    -q2(2)*q1(3);
            qn = [qnew_0;qnew_1;qnew_2;qnew_3];
        end

        function qmat = hamilton_product_as_matrix(q)
            qmat = [...
                q(1), -q(2), -q(3), -q(4)
                q(2),  q(1), -q(4), q(3)
                q(3), q(4), q(1), -q(2)
                q(4), -q(3), q(2), q(1)];
        end

        function data = fuse_data(data_sources, noises)
            % Mainly for building the u vector if some quantity
            % is captured by more than one sensor (ex: redundancy)
            data = ivwmean(data_sources, noises);
        end

        function rot_mat = rotmat(q)
            % Returns a quaternion's equivalent rotation matrix
            q0 = q(1);
            q1 = q(2);
            q2 = q(3);
            q3 = q(4);

            rot_mat = [...
                q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)
                2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1)
                2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3];
        end

    end
end