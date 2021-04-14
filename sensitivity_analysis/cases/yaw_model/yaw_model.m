function Y = yaw_model(X,P)
% this is the model for the induced velocities in yaw as proposed by
% Schepers in his PhD thesis (Appendix B)

% processing the parameters
switch nargin
    case 1
        phi_y = 0;
        phi_r = 0;
    case 2
        phi_y = P(1);
        phi_r = P(2);
    otherwise
        error('Number of input arguments not accepted!');
end

%
% n_rel = 5;
% r_rel = linspace(0,1,n_rel);

r_rel = 0.5;

% X contains 20 parameters
f0_A1 = X(:,1);
f1_A1 = X(:,2);
f2_A1 = X(:,3);
f3_A1 = X(:,4);
f4_A1 = X(:,5);
f0_A2 = X(:,6);
f1_A2 = X(:,7);
f2_A2 = X(:,8);
f3_A2 = X(:,9);
f4_A2 = X(:,10);
f0_psi1 = X(:,11);
f1_psi1 = X(:,12);
f2_psi1 = X(:,13);
f3_psi1 = X(:,14);
f4_psi1 = X(:,15);
f0_psi2 = X(:,16);
f1_psi2 = X(:,17);
f2_psi2 = X(:,18);
f3_psi2 = X(:,19);
f4_psi2 = X(:,20);

if (abs(phi_y)>=15)
    A1 = Ampl1(r_rel,phi_y);
    A2 = Ampl2(r_rel,phi_y);
else
    A1 = (abs(phi_y)/15) * Ampl1(r_rel,phi_y);
    A2 = (abs(phi_y)/15) * Ampl2(r_rel,phi_y);
end

if (r_rel <= 0.15)
    psi1 = 270;
elseif (r_rel > 0.15 && r_rel < 0.29)
    psi1 = 270 + ((r_rel - 0.15)/0.14)*(Shift1(0.29,phi_y) - 270);
elseif (r_rel >= 0.29 && r_rel <= 0.71)
    psi1 = Shift1(r_rel,phi_y);
elseif (r_rel > 0.71 && r_rel < 0.85)
    psi1 = Shift1(0.71,phi_y) + ((r_rel-0.71)/0.14)*(90 - Shift1(0.71,phi_y));
elseif (r_rel >= 0.85)
    psi1 = 90;
else
    error('wrong value for r_rel');
end


if (abs(phi_y)>=0 && abs(phi_y)<=90)
    psi2 = Shift2(r_rel,phi_y);
else
    error('this value of phi_y is not supported by the yaw model');
end
    

% compute the response value
Y(:,1) = 1 - A1*cos(phi_r - psi1) - A2*cos(2*phi_r - psi2);


    function A1 = Ampl1(r_rel,phi_y)
        A1 = f0_A1 + f1_A1.*r_rel + f2_A1*r_rel.^2 + f3_A1*sind(abs(phi_y)) + ...
            f4_A1*sind(phi_y)^2;
    end

    function A2 = Ampl2(r_rel,phi_y)
        A2 = f0_A2 + f1_A2.*r_rel + f2_A2*r_rel.^2 + f3_A2*sind(abs(phi_y)) + ...
            f4_A2*sind(phi_y)^2;
    end
    function Psi1 = Shift1(r_rel,phi_y)
        Psi1 = f0_psi1 + f1_psi1.*r_rel + f2_psi1*r_rel.^2 + f3_psi1*sind(abs(phi_y)) + ...
            f4_psi1*sind(phi_y)^2;
    end
    function Psi2 = Shift2(r_rel,phi_y)
        Psi2 = f0_psi2 + f1_psi2.*r_rel + f2_psi2*r_rel.^2 + f3_psi2*sind(abs(phi_y)) + ...
            f4_psi2*sind(phi_y)^2;
    end
end