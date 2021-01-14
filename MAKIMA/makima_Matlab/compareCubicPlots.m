%% |compareCubicPlots|
function compareCubicPlots(x,v,xq,showMakima,legendLocation)
% Plot 'pchip', 'spline', Akima, and 'makima' interpolation results
    vqp = pchip(x,v,xq);      % same as vq = interp1(x,v,xq,'pchip')
    vqs = spline(x,v,xq);     % same as vq = interp1(x,v,xq,'spline')
    vqa = akima(x,v,xq);
    if showMakima
        vqm = makima(x,v,xq); % same as vq = interp1(x,v,xq,'makima')
    end

    plot(x,v,'ko','LineWidth',2,'MarkerSize',10,'DisplayName','Input data')
    hold on
    plot(xq,vqp,'LineWidth',4,'DisplayName','''pchip''')
    plot(xq,vqs,'LineWidth',2,'DisplayName','''spline''')
    plot(xq,vqa,'LineWidth',2,'DisplayName','Akima''s formula')
    if showMakima
        plot(xq,vqm,'--','Color',[0 3/4 0],'LineWidth',2, ...
            'DisplayName','''makima''')
    end
    hold off
    xticks(x(1)-1:x(end)+1)
    legend('Location',legendLocation)
    title('Cubic Hermite interpolation')
end


