function sng_chistplot(Eh,Bh,minC);
%% computes the overlap en plot both histograms and the overlap in one image

    figure
    for i = 1:size(Eh,1)     
        OL(i,:)=min(Eh(i,:),Bh(i,:));
        subplot(size(Eh,1),1,i);
        hold on;
            %set(gca,'Xlim',[0,])
            bar(Eh(i,:),'Facecolor',[0.4,0.4,0.8],'Edgecolor','none','BarWidth',1);set(gca,'XLim',[0 length(Eh)])
            bar(Bh(i,:),'Facecolor',[0.4,0.8,0.4],'Edgecolor','none','BarWidth',1);set(gca,'XLim',[0 length(Eh)])
            bar(OL(i,:),'Facecolor',[0.2,0.6,0.6],'Edgecolor','none','BarWidth',1);set(gca,'XLim',[0 length(Eh)])
            line([minC(i) minC(i)],[get(gca,'Ylim')],'color',[0 0 0]);
            %annotation(figure1,'line',[0.5660 0.56607],[0.856142857142857 0.707142857142857]);
        hold off;
    end    
end