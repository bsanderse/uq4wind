function [Cv,xp,yp] = get_chord(x,Rt,ig,color)

xcp = [0 1.1439935563517057e+01   2.3345657075283732e+01   3.5488860834847088e+01   3.8383411004288320e+01 Rt];

ycp = x(1:end);

xp = [0
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

yp = interp1(xcp,ycp,xp,'pchip');

a = sign(diff(yp)); a(diff(a)==0)=[]; sign_ch_yp = length(a);

if ig==-1
    figure(2)
    set(gca,'Fontsize',16);
    plot(xcp,ycp,'o','Color',color,'LineWidth',2)
    hold on
    plot(xp,yp,'-','Color',color,'LineWidth',2)
    hold on
    grid on
    axis([0 65 0 7])
    xlabel('r [m]');
    ylabel('c [m]');
end

% if sign_ch_yp > 2
%     Cv = sign_ch_yp;
% else
%     Cv = 0;
% end

Cv = 0;