function figlabels(xlabel,ylabel,zlabel, figtitle,size)

    xlab = get(gca,'XLabel'); set(xlab,'String',xlabel,'FontSize',size,'FontName','Arial')
    ylab = get(gca,'YLabel'); set(ylab,'String',ylabel,'FontSize',size,'FontName','Arial')
    zlab = get(gca,'ZLabel'); set(zlab,'String',zlabel,'FontSize',size,'FontName','Arial')
    title(figtitle,'FontSize',size,'FontWeight','Bold')

end