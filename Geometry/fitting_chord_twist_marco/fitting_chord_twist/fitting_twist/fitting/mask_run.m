function [J] = mask_run(x)

color = 'b';

ig = 1;

Rt = 38.8; % tip radius

[Cv,xp,yp] = get_twist(x,Rt,ig,color);

if Cv==0
    
    RNodes = [0
        2
        4
        6
        8
        10
        12
        14
        16
        18
        20
        22
        24
        26
        28
        30
        32
        34
        36
        37
        38
        38.4
        38.8];
    
    ref = [0
        5.37
        6.69
        7.9
        9.11
        10.19
        9.39
        7.16
        5.45
        4.34
        3.5
        2.86
        2.31
        1.77
        1.28
        0.9
        0.55
        0.23
        0.03
        0.02
        0.93
        2.32
        6.13];
    
    if ig==-1
        figure(1)
        box on
        set(gca,'Fontsize',16);
        plot(xp,yp,'-b','LineWidth',2)
        hold on
        plot(RNodes,ref,'-r','LineWidth',2)
        hold off
        grid on
        axis([0 65 0 25])
        xlabel('r [m]');
        ylabel('\theta [deg]');
    end
    
    n = length(ref(1:end,1));
    J = sqrt((sum((yp(1:end,1)-ref(1:end,1)).^2))/n); % Root-mean-square deviation
    
else
    
    J = 1000000+(Cv);
    
end





