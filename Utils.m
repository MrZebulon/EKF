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

        function mat = cross_product_matrix(vec)

            
            mat = [...
                0, -vec(3), -vec(2)
                vec(3), 0, -vec(1)
                vec(2), vec(1), 0];

        end
    end
end