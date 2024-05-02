% Funções para matéria de robótica
classdef frb
   methods(Static)
       % metodos para calcular o rotacional simbolico de uma matriz
       %[x y z] [n o a] 
       function X = srotx(theta)
           X = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
       end
       function X = sroty(theta)
           X = [cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
       end
       function X = srotz(theta)
           X = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
       end
       
       function X = Rot(theta, x)
            if (isnumeric(theta))
                switch x
                    case 'x'
                        X = [rotx(theta) [0 0 0].'; zeros(1,3) 1];
                    case 'y'
                        X = [roty(theta) [0 0 0].'; zeros(1,3) 1];
                    case 'z'
                        X = [rotz(theta) [0 0 0].'; zeros(1,3) 1];
                end
            else
                switch x
                    case 'x'
                        X = [frb.srotx(theta) [0 0 0]'; zeros(1,3) 1];
                    case 'y'
                        X = [frb.sroty(theta) [0 0 0]'; zeros(1,3) 1];
                    case 'z'
                        X = [frb.srotz(theta) [0 0 0]'; zeros(1,3) 1];
                end
                
            end
       end
       
       function X = Trans(P)
            X = [eye(3) P.'; zeros(1,3) 1];
        end
                     
       function X = T_cart(P)
            X = rb.Trans(P);
        end
        
       function X = T_cilin(r,a,l)
            X = frb.Trans([0 0 l])*frb.Rot(a,'z')*frb.Trans([r 0 0]);
        end

       function X = T_esfe(r,b,y)
            X = frb.Rot(y,'z')*frb.Rot(b,'y')*frb.Trans([0 0 r]);
        end
        
       function [px py pz] = pos_cart(M)
            arguments M (4,4)
            end
            px = M(1,4);
            py = M(2,4);
            pz = M(3,4);
        end
        
       function [r alpha l] = pos_cil(M)
            arguments M (4,4)
            end
            alpha = atan2d(M(2,4),M(1,4));
            r = M(2,4)/sind(alpha);
            l = M(3,4);
        end
        
       function [r beta gamma] = pos_esf(M)
            arguments M (4,4)
            end
            gamma = atan2d(M(2,4),M(1,4));
            rsinB = M(2,4)/sind(gamma);
            beta = atan2d(rsinB,M(3,4));
            r = M(3,4)/cosd(beta);
        end
        
       function [n o a] = t2rag(M)
            arguments M (4,4)
            end
            % movimentos fixos 
            % XYZ
            % phi a
            % theta o
            % psi n
            a = atan2d(M(2,1),M(1,1));
            o = atan2d(-M(3,1),(M(1,1)*cosd(a)+ M(2,1)*sind(a)));
            n = atan2d((-M(2,3)*cosd(a)+ M(1,3)*sind(a)),(M(2,2)*cosd(a)-M(1,2)*sind(a)));
 
       end
             
       function X = rag2t(n,o,a)
           X = frb.Rot(a,'z')*frb.Rot(o,'y')*frb.Rot(n,'x');
       end
       
       function X = euler2t(phi,theta,psi)
           X = frb.Rot(phi,'z')*frb.Rot(theta,'y')*frb.Rot(psi,'z');
       end
       
       function [phi theta psi somaphipsi subphipsi] = t2euler(R, sequence)
            % Euler movimentos correntes
            % rad2deg()
            % tform2eul(M,'ZYZ')
            % Extract the rotation matrix elements
            r11 = R(1,1);
            r12 = R(1,2);
            r13 = R(1,3);
            r21 = R(2,1);
            r22 = R(2,2);
            r23 = R(2,3);
            r31 = R(3,1);
            r32 = R(3,2);
            r33 = R(3,3);
            subphipsi = nan;
            somaphipsi = nan;
            phi = nan;
            psi = nan;
            
            % Determine the sequence of rotations and calculate Euler angles
            switch sequence
                case 'ZYZ' 
                    if r33 == 1
                        theta = 0;                       
                        somaphipsi = atan2d(r21,r11);
                    else
                        theta = atan2d(sqrt(1-r33^2),r33);
                        phi = atan2d(r23,r13);
                        psi = atan2d(r32,-r31);
                    end
                case 'XYZ' 
                    if r13 == 1
                        theta = 90;
                        somaphipsi = atan2d(r32,r22);
                    else
                        theta = atan2d(r13,sqrt(1-r13^2));
                        phi = atan2d(-r23,r33);
                        psi = atan2d(-r12,r11);
                    end
                case 'XYX' 
                    if r11 == 1
                        theta = 0;
                        somaphipsi = atan2d(r32,r22);
                    else
                        theta = atan2d(sqrt(1-r11^2),r11);
                        phi = atan2d(r21,-r31);
                        psi = atan2d(r12,r13);
                    end
                case 'XZX' 
                    if r11 == 1
                        theta = 0;
                        somaphipsi = atan2d(r32,r22);
                    else
                        theta = atan2d(sqrt(1-r11^2),r11);
                        phi = atan2d(r31,r21);
                        psi = atan2d(r13,-r12);
                    end
                case 'YXY' 
                    if r22 == 1
                        theta = 0;
                        somaphipsi = atan2d(r13,r33);
                    else
                        theta = atan2d(sqrt(1-r22^2),r22);
                        phi = atan2d(r12,r32);
                        psi = atan2d(r21,-r23);
                    end
                case 'YZY' 
                    if r22 == 1
                        theta = 0;
                        somaphipsi = atan2d(r13,r33);
                    else
                        theta = atan2d(sqrt(1-r22^2),r22);
                        phi = atan2d(r32,-r12);
                        psi = atan2d(r23,r21);
                    end
                case 'ZXZ' 
                    if r33 == 1
                        theta = 0;
                        somaphipsi = atan2d(r21,r11);
                    else
                        theta = atan2d(sqrt(1-r33^2),r33);
                        phi = atan2d(r13,-r23);
                        psi = atan2d(r31,r32);
                    end
                    % Cardaniano
                case 'XZY' 
                    if r12 == -1
                        theta = 90;
                        subphipsi = atan2d(r31,r33);
                    else
                        theta = atan2d(-r12,sqrt(1-r12^2));
                        phi = atan2d(r32,r22);
                        psi = atan2d(r13,r11);
                    end
                case 'YZX' 
                    if r21 == 1
                        theta = 90;
                        somaphipsi = atan2d(r13,r33);
                    else
                        theta = atan2d(r21,sqrt(1-r21^2));
                        phi = atan2d(-r31,r11);
                        psi = atan2d(-r23,r22);
                    end
                case 'YXZ' 
                    if r23 == -1
                        theta = 90;
                        subphipsi = atan2d(r12,r11);
                    else
                        theta = atan2d(-r23,sqrt(1-r23^2));
                        phi = atan2d(r13,r33);
                        psi = atan2d(r21,r22);
                    end
                case 'ZXY' 
                    if r32 == 1
                        theta = 90;
                        somaphipsi = atan2d(r21,r11);
                    else
                        theta = atan2d(r32,sqrt(1-r32^2));
                        phi = atan2d(-r12,r22);
                        psi = atan2d(-r31,r33);
                    end
                case 'ZYX' 
                    if r31 == -1
                        theta = 90;
                        subphipsi = atan2d(r23,r13);
                    else
                        theta = atan2d(-r31,sqrt(1-r31^2));
                        phi = atan2d(r21,r11);
                        psi = atan2d(r32,r33);
                    end
                otherwise
                    error('Unsupported rotation sequence.');
            end

       end
       
       function X = t2pc(M06, d6)
           % pc = d06 - d6*R06*
           X = M06(1:3,4) - d6*M06(1:3,1:3)*[0 0 1].';
       end
       
       function X = t2ori(M06,R03)
           %retorna R36
           % R36 = inv(r03)*R06
           X = inv(R03)*M06(1:3,1:3);
       end
       function X = denavit(a,alfa,d,theta)         
          X =  frb.Rot(theta,'z')*frb.Trans([a 0 0])*frb.Trans([0 0 d])*frb.Rot(alfa,'x');        
       end
       
       function X = jacobiano(tipoJuntas, DH)
            T = eye(4);
            dh = eye(4);
            n = length(DH)/4;
            J = zeros(6,n);
%             if ~isnumeric(DH)                
%                 dh = sym(dh);
%                 T = sym(T);
%                 J = sym(J);
           % end           
            a = 1;
            for i = 1:n
                dh(:,:,i) = DH(1:4,a:a+3);
                T = T * dh(:,:,i);
                a = a + 4;
            end
            
            
            z0 = [0 0 1]'

            
            if strcmp(tipoJuntas(1),"r")
                O0 = [0 0 0]'
                On = T(1:3,4);
                fprintf("O%i\n",n)
                disp(On)
                J1 = [cross(z0,On-O0); z0];
            else
                J1 = [z0; zeros(3,1)];
            end
            J(:,1) = J1;
            
            for i = 1:n-1
                A = eye(4);
                for a = 1:i
                    A = A* dh(:,:,a);
                end    
                zi_1 = A(1:3,1:3)*[0 0 1]';
                fprintf("Z%i\n",i)
                disp(zi_1)
                
                if strcmp(tipoJuntas(i+1),"r")
                    Oi_1 = A(1:3,4);
                    fprintf("O%i\n",i)
                    disp(Oi_1)
                    Ji = [cross(zi_1,On-Oi_1); zi_1];
                else
                    Ji = [zi_1; zeros(3,1)];
                end
                J(:,i+1) = Ji;
            end
            X = J;
       end
       
       function [c0 c1 c2 c3] = coef3(theta1, theta2, tf)
           B = [theta1, 0, theta2, 0]';
           A = [1 0 0 0; 0 1 0 0; 1 tf tf^2 tf^3; 0 1 2*tf 3*tf^2];
           X = A\B;
           c0 = X(1);
           c1 = 0;
           c2 = X(3);
           c3 = X(4);
           syms T
           pos = vpa(c0  +c2*T^2 + c3*T^3,4)
           vel = vpa( 2*c2*T + 3*c3*T^2,4)
           acel = vpa(2*c2 + 6*c3*T,4)
       end 
       
       
       function [c0 c1 c2 c3 c4 c5] = coef5(theta1, theta2, tf, acel)
           c0 = theta1;
           c1 = 0;
           c2 = acel/2;
           B = [theta2-c0-c1*tf-c2*tf^2,-c1-2*c2*tf, -acel-2*c2]';
           A = [tf^3 tf^4 tf^5; 3*tf^2 4*tf^3 5*tf^4; 6*tf 12*tf^2 20*tf^3];
           X = A\B;
           c3 = X(1);
           c4 = X(2);
           c5 = X(3);
           syms T
           pos = vpa(c0 +c2*T^2 + c3*T^3 + c4*T^4 + c5*T^5,4)
           vel = vpa( 2*c2*T + 3*c3*T^2 + 4*c4*T^3 + 5*c5*T^4,4)
           acel = vpa(2*c2 + 6*c3*T + 12*c4*T^2 + 20*c5*T^3,4)
           
       end
       
       function X = matrixx(t1f, t2f, t3f)
  
            L1 = [1, zeros(1,13)];
            L2 = [0, 1, zeros(1,12)];
            L3 = [0, 0, 2, zeros(1,11)];
            L4 = [1, t1f, t1f^2, t1f^3, t1f^4, zeros(1,9)];
            L5 = [zeros(1,5), 1, zeros(1,8)];
            L6 = [0, 1, 2*t1f, 3*t1f^2, 4*t1f^3, 0, -1, zeros(1,7)];
            L7 = [0, 0, 2, 6*t1f, 12*t1f^2, 0, 0, -2, zeros(1,6)];
            L8 = [zeros(1,5), 1, t2f, t2f^2, t2f^3, zeros(1,5)];
            L9 = [zeros(1,9), 1, zeros(1,4)];
            L10 = [zeros(1,6), 1, 2*t2f, 3*t2f^2, 0, -1,0, 0, 0];
            L11 = [zeros(1,7), 2, 6*t2f, 0, 0, -2, 0, 0];
            L12 = [zeros(1,9), 1, t3f, t3f^2, t3f^3, t3f^4];
            L13 = [zeros(1,10), 1, 2*t3f, 3*t3f^2, 4*t3f^3];
            L14 = [zeros(1,11), 2, 6*t3f, 12*t3f^2];

            X = [L1; L2; L3; L4; L5; L6; L7; L8; L9; L10; L11; L12; L13; L14];
         
       end
       
   end
end