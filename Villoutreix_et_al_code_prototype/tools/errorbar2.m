function [hh] = errorbar2(x,y,l,col),
% ERRORBAR2 generates a plot with colored error bar - x contains the 
% abscissa - y contains the ordinates - l contains the error - col contains
% the color of the points

if ((length(x) == length(y))&(length(x) == length(l))),
    for k = 1:length(x),
        hold on,
        plot(x(k),y(k),'o','MarkerFaceColor',col,'MarkerSize',10,'MarkerEdgeColor',col);
        hold off,
        tee = 0.075;
        xl = x(k) - tee;
        xr = x(k) + tee;
        if l(k)>0,
        ytop = y(k) + l(k);
        ybot = y(k) - l(k);
        step = (ytop-ybot)/100;
        end
        hold on,
        plot(xl:0.01:xr,ytop*ones(length(xl:0.01:xr),1),'LineWidth',4,'Color',col);
        hold on,
        plot(xl:0.01:xr,ybot*ones(length(xl:0.01:xr),1),'LineWidth',4,'Color',col);
        hold on, 
        plot(x(k)*ones(length(ybot:step:ytop),1),ybot:step:ytop,'LineWidth',4,'Color',col);
    end
end
end