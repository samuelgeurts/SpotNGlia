function sng_figcm(fsx,fsy,DPI)
%gives image real size in centimeters use for export with the export_fig function
%because the dpi is fixed to 72, the screen size differs per screen
%Give the real dpi value to create a copy of the figure such that it is
%displayed with real size in centimeters

%Example: sng_figcm(10,8,113.6)

%{
fsx = obj.fsxy(1)
fsy = obj.fsxy(2)
%}

    h1 = gcf;
    set(h1,'PaperUnits','centimeters','Color',[1 1 1]);
    set(h1,'Units','centimeters');
    pos = get(h1,'Position');
    pos(3) = fsx;
    pos(4) = fsy;
    set(h1,'Units','centimeters','Position',pos);

    if exist('DPI','var')
        ScaledFigure.calibrateDisplay(DPI);
        ScaledFigure(gcf,'copy');
        set(gcf,'Units','Centimeters');
        set(gcf,'Position',(get(gcf,'Position') + [fsx 0 0 0]));

    end
    
    
    
end
