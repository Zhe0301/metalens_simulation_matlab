function plot_phase(L, varargin)
    p = inputParser;            
    addParameter(p,'show',false);      
    addParameter(p,'save_path','./'); 
    parse(p,varargin{:}); 
    show = p.Results.show;
    save_path = p.Results.save_path;
    
    phase_2pi = mod(L.phase, 2 * pi);
    if show
        f = figure;
    else
        f = figure(Visible="off");
    end
    f.Position(3:4) = [1024 512];
    subplot(1, 2, 1);
    
    imagesc(L.Grid.d2_x(1,:), L.Grid.d2_y(:,1), L.phase);
    ax=gca;
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    ax.DataAspectRatio = [1,1,1];
    title('Phase Distribution',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{x}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");
    c = colorbar;
    c.Label.String = 'Phase (rad)';
    c.Label.FontSize = 16;
    c.Label.FontName = "Times New Roman";

    subplot(1, 2, 2);    
    imagesc(L.Grid.d2_x(1,:), L.Grid.d2_y(:,1), phase_2pi);
    ax=gca;
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    ax.DataAspectRatio = [1,1,1];
    title('Phase Distribution 0 to 2\pi',FontSize=16,FontName="Times New Roman",Interpreter="tex");
    xlabel('{\it{x}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");
    c = colorbar;
    c.Label.String = 'Phase (rad)';
    c.Label.FontSize = 16;
    c.Label.FontName = "Times New Roman";


    print(f,fullfile(save_path, 'Phase.png'),'-dpng','-r300')
    if show
        waitfor(f)
    else
        close(f);
    end
end