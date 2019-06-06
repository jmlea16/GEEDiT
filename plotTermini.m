function plotTermini(Results)
% figure('Color',[1,1,1])
% propertyeditor('on')
% plotbrowser('on')

 
 figure('units','normalized','outerposition',[0 0 1 1],'Color',[1,1,1])
 figure(1)
% if strcmp(Results.Method,'Curvilinear Box Method')||strcmp(Results.Method,'Variable Box Method')
subplot(2,2,1)
%  end
plot(Results.Centreline(:,1),Results.Centreline(:,2),'--','DisplayName','Centreline','Color','black')
hold on
colour=colormap(jet(length(Results.Date(:,1))));
for n=1:length(Results.TerminusGeometry(:,1))
    if ~isempty(Results.TerminusGeometry{n,1})
    plot(Results.TerminusGeometry{n,1}(1,:),Results.TerminusGeometry{n,1}(2,:),...
        'DisplayName',strcat(num2str(Results.Date(n,3)),'/',num2str(Results.Date(n,2)),'/',...
        num2str(Results.Date(n,1))),'Color',colour(n,:))
    text(Results.TerminusGeometry{n,1}(1,end-1),Results.TerminusGeometry{n,1}(2,end-1),...
        cellstr(num2str(round(Results.Distance(n,1),0))),'Color',colour(n,:))
    end
end
% legend('show')
% set(legend,'Location','northeastoutside');
title('Glacier termini and their distances from the most recent observation')
xlabel(Results.Method)
 axis equal
hold off
% end

if strcmp(Results.Method,'Curvilinear Box Method')||strcmp(Results.Method,'Variable Box Method')
 figure(1)
subplot(2,2,2)
plot(Results.Centreline(:,1),Results.Centreline(:,2),'--','DisplayName','Centreline','Color','black')
hold on
for n=1:length(Results.TerminusGeometry(:,1))
    if ~isempty(Results.TerminusGeometry{n,1})
    plot(Results.BoxGeometry{n,1}(:,1),Results.BoxGeometry{n,1}(:,2),...
        'DisplayName',strcat(num2str(Results.Date(n,3)),'/',num2str(Results.Date(n,2)),'/',...
        num2str(Results.Date(n,1))),'Color',colour(n,:))
    text(Results.TerminusGeometry{n,1}(1,end-1),Results.TerminusGeometry{n,1}(2,end-1),...
        cellstr(num2str(round(Results.Distance(n,1),0))),'Color',colour(n,:))
    end
end
title('Boxes used to calculate terminus change')
xlabel(Results.Method)
% legend('show')
 axis equal
% set(legend,'Location','northeastoutside');
hold off
end

 figure(1)
subplot(2,2,3:4)
plot(Results.Date(:,4),Results.Distance(:,1),'-x','DisplayName','Observation');
% legend('show')
% set(legend,'Location','northeastoutside');
title('Terminus position relative to most recent observation (m)')
datetick()
xlabel('Year')
ylabel('Terminus position (m')
grid on
reference=refline(0,0)
reference.Color='black';
hold off



end