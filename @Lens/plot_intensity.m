function plot_intensity(L, varargin)
    p = inputParser;            
    addParameter(p,'show',false);      
    addParameter(p,'save_path','./'); 
    parse(p,varargin{:}); 
    show = p.Results.show;
    save_path = p.Results.save_path;
    if show
        f = figure;
    else
        f = figure(Visible="off");
    end
    f.Position(3:4) = [512 512];
    imagesc(L.Grid.d2_x(1,:), L.Grid.d2_y(:,1), abs(L.complex_amplitude_t));
    ax = gca;
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    ax.DataAspectRatio = [1,1,1];
    title('Intensity',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{x}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");
    c = colorbar;
    c.Label.String = 'Intensity';
    c.Label.FontSize = 16;
    c.Label.FontName = "Times New Roman";
    exportgraphics(f,fullfile(save_path, 'Intensity.png'),'Resolution',300)
    if show
        waitfor(f)
    else
        close(f);
    end

end