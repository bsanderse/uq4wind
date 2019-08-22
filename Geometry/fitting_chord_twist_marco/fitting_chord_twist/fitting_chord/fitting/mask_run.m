function [J] = mask_run(x)

color = 'b';

ig = 1;

Rt = 38.8; % tip radius

[Cv,xp,yp] = get_chord(x,Rt,ig,color);

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
    
    ref = [2.42
        2.48
        2.65
        2.81
        2.98
        3.14
        3.17
        2.99
        2.79
        2.58
        2.38
        2.21
        2.06
        1.92
        1.8
        1.68
        1.55
        1.41
        1.18
        0.98
        0.62
        0.48
        0.07];
    
    if ig==-1
        figure(1)
        box on
        set(gca,'Fontsize',16);
        plot(xp,yp,'-b','LineWidth',2)
        hold on
        plot(RNodes,ref,'-r','LineWidth',2)
        hold off
        grid on
        axis([0 65 0 5])
        xlabel('r [m]');
        ylabel('c [m]');
    end
    
    n = length(ref(1:end,1));
    J = sqrt((sum((yp(1:end,1)-ref(1:end,1)).^2))/n); % Root-mean-square deviation
    
else
    
    J = 1000000+(Cv);
    
end





