%% settings figure

fsx = 8.5
fsy = 8.5
h1 = figure('PaperUnits','centimeters','Color',[1 1 1]);
sng_figcm(fsx,fsy)

%% plot axes

hold on
line([-10,4],[0,0],'Color',[0 0 0],'LineWidth',0.5)
line([0,0],[-4,10],'Color',[0 0 0],'LineWidth',0.5)
%s1 = scatter(WBR(:,1),WBR(:,2),'MarkerEdgeColor',[0 0.7 0 ],'SizeData',10,'LineWidth',0.5,'Marker','o')
%s2 = scatter(WGR(:,1),WGR(:,2),'MarkerEdgeColor',[0 0 0.7 ],'SizeData',10,'LineWidth',0.5,'Marker','o')
%l=legend([s1,s2],'green to red','blue to red','Location','northwest')

%% settings axis

xlabel('x-displacement [pix]','FontSize',8,'FontName','arial');
ylabel('y-displacement [pix]','FontSize',8,'FontName','arial');

set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
set(gca,'XLim',[-10,4])
set(gca,'YLim',[-4,10])
set(gca,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])

%% show real size on screen

sng_figcm(fsx,fsy,113.6)
set(0, 'currentfigure', h1)

%% export
export_fig(h1 ,['/Users/samuelgeurts/Desktop/','fig',num2str(1)], '-png', '-r600', '-nocrop');