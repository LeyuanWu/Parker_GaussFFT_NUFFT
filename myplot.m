function plot1=myplot(ax,X,Ys,fs,MS,dispNames,begInd)
%%%% Red,Blue,Cyan,Magenta,Green,Orange,Pruple,Yellow
Colors=[255 0 0;0,0,255;0 255 255;255 0 255;
    0 255 0;255 97 0;160 32 240;0 139 0;255 255 0]/255;
Markers={'o','^','+','>','x','<','v','diamond','square','pentagram'};
lineStyles={'none','none','none',':',':',':',':',':',':',':'};
%%%% Mapping
hold(ax,'all');grid(ax,'on');
set(ax,'Fontsize',fs);
plot1=plot(ax,X,Ys);
for i=1:1:size(Ys,2)
    set(plot1(i),'Marker',Markers{i+begInd},'Color',Colors(i+begInd,:),...
        'LineStyle',lineStyles{i+begInd},'MarkerSize',MS(i+begInd),...
        'MarkerFaceColor','none',...
        'DisplayName',dispNames{i});
end
legend('Location','best')
end

