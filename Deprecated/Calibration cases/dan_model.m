function Y = dan_model(Q,P)
% Q: uncertain parameter, size p*2
% P: contains independent variables, P = [beta c V x]
% ind: range of each model
Y = zeros(4,1);

Y(1) = 0.5*1.225.*P(1,2).*(P(1,3).^2).*( ...
    Q(1).*sin(P(1,1)) + Q(2).*cos(P(1,1)) );

Y(2) = 0.5*1.225.*P(2,2).*(P(2,3).^2).*( ...
    Q(3).*sin(P(2,1)) + Q(4).*cos(P(2,1)) );

Y(3) = 0.5*1.225.*P(3,2).*(P(3,3).^2).*( ...
    Q(5).*sin(P(3,1)) + Q(6).*cos(P(3,1)) );

Y(4) = 0.5*1.225.*P(4,2).*(P(4,3).^2).*( ...
    Q(7).*sin(P(4,1)) + Q(8).*cos(P(4,1)) );

