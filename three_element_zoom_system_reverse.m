function [e_5, e_6, e_yz] = three_element_zoom_system_reverse(S, L1, L2, L3, G, d_lens, d_12, d_23, d_bfl, efl,varargin)
    % 3片式变焦系统, 超透镜位于每片的后表面
    % S: 光源对象
    % L1, L2, L3: 三片镜片对象
    % G: 网络对象
    % d_lens, d_12, d_23, d_bfl: 镜片厚度 1，2片间距 2，3片间距 后截距
    % efl 有效焦距，用于文件命名
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 可选参数
    %··refractive_index 镜片折射率
    % save_path 文件储存目录
    % magnification 图像放大倍率，为0或小于0时不放大,算y-z平面光场（sampling_point采样点数,interval焦点前后距离，单位mm，二者之一为0或小于0时不计算）
    % sampling_point x-z截面计算的采样点数，为0时不计算
    % interval x-z截面采样时焦点前后的距离，为0时不计算，单位mm
    p = inputParser;    
    addParameter(p,'refractive_index',1); 
    addParameter(p,'save_path','./');
    addParameter(p,'sampling_point',0);
    addParameter(p,'interval',0);
    addParameter(p,'method','AS');
    addParameter(p,'show',false);
    addParameter(p,'gpu_acceleration',false);
    parse(p,varargin{:});       
    refractive_index = p.Results.refractive_index;
    save_path = p.Results.save_path;
    method = p.Results.method;
    sampling_point = p.Results.sampling_point;
    interval = p.Results.interval;
    show = p.Results.show;
    gpu_acceleration = p.Results.gpu_acceleration;

    if gpu_acceleration
        disp("Using GPU to accelerate computing");
    end
    t0 = tic;
    
    %% 参数初始化
    p_b = PropOperator(G, S.wavelength_vacuum, d_lens, 'refractive_index', refractive_index, 'method', method, 'gpu_acceleration', gpu_acceleration);
    p_12 = PropOperator(G, S.wavelength_vacuum, d_12, 'method', method, 'gpu_acceleration', gpu_acceleration);
    p_23 = PropOperator(G, S.wavelength_vacuum, d_23, 'method', method, 'gpu_acceleration', gpu_acceleration);
    p_bfl = PropOperator(G, S.wavelength_vacuum, d_bfl, 'method', method, 'gpu_acceleration', gpu_acceleration);
    
    elapsed_time = toc(t0);
    fprintf("Parameter initialization is completed. Elapsed time: %.2f s\n", elapsed_time);
    
    %% 光传播计算
    e_1 = p_b.prop(S.complex_amplitude);
    e_2 = p_12.prop(e_1 .* L1.complex_amplitude_t);
    e_3 = p_b.prop(e_2);
    e_4 = p_23.prop(e_3 .* L2.complex_amplitude_t);
    e_5 = p_b.prop(e_4);
    e_5 = e_5 .* L3.complex_amplitude_t;
    e_6 = p_bfl.prop(e_5);
   
    %% 像面绘图
    % 计算光强和中间值
    I = abs(e_6).^2;
    mid_index_0 = round(length(G.axis) / 2);
    if show
        f = figure;
    else
        f = figure(Visible="off");
    end
    f.Position = [0 0 1280 1024];
    subplot(2, 2, 1);
   
    imagesc(G.d2_x(1,:), G.d2_y(:,1), I);
    ax=gca;
    ax.DataAspectRatio = [1,1,1];
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    
    c = colorbar;
    c.Label.String = 'Intensity';
    c.Label.FontSize = 16;
    c.Label.FontName = "Times New Roman";
    title('Image',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{x}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");
    
    subplot(2, 2, 2);
    plot(G.axis, I(mid_index_0, :)); 
    ax=gca;
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    title('{\it{x}}-axis Cross Section',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{x}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('Intensity',FontSize=16,FontName="Times New Roman");
    grid on;
    
    subplot(2, 2, 3);
    imagesc(G.d2_x(1,:), G.d2_y(:,1), log10(I));
    ax=gca;
    ax.DataAspectRatio = [1,1,1];
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    ax.DataAspectRatio = [1,1,1];
    c = colorbar;
    c.Label.String = 'log10(Intensity)';
    c.Label.FontSize = 16;
    c.Label.FontName = "Times New Roman";
    title('Image',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{x}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");

    subplot(2, 2, 4);
    plot(G.axis, log10(I(mid_index_0, :)));
    ax=gca;
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    title('{\it{y}}-axis Cross Section'); 
    xlabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");
    ylabel('log10(Intensity)',FontSize=16,FontName="Times New Roman");
    grid on;
    exportgraphics(f,fullfile(save_path,sprintf('image_f_%.1f.png', efl)),'Resolution',300)

    if show
        waitfor(f)
    else
        close(f);
    end

    %% 计算FWHM和最大光强位置
    [center_x, center_y, FWHM] = calculate_FWHM(I, G);
    FWHM = FWHM*2;
    fprintf("centriod = (%.4f, %.4f) mm\n",center_x,center_y);
    fprintf("FWHM = %.2f μm\n",FWHM*1000);
    %% 圈入能量计算
    [radius, energy_ratio] = calculate_enclosed_energy_ratio(I, G, ...
        center_x, center_y,'max_radius',12e-3);
    
    F = griddedInterpolant(radius,energy_ratio,'pchip'); % 三次样条插值
    radius_interp = linspace(min(radius), max(radius), 3 * length(radius)); % 插值2~3倍
    energy_ratio_interp = F(radius_interp);
    if show
        f = figure;
    else
        f = figure(Visible="off");
    end
    f.Position = [0 0 640 512];
    plot(radius_interp*1e3, energy_ratio_interp,'Color',"#4169E1",'LineWidth',2);  % #4169E1皇家蓝
    ax = gca;
    ax.FontSize = 16;
    ax.FontName = "Times New Roman";
    ax.XLim = [0 10];
    ax.YLim = [0 1];
    grid on
    hold on;
    l1 = xline(FWHM/2*1e3,'--', 'Color',"#FF0000",'LineWidth',2); % #FF0000红色
    l2 = xline(2* FWHM/2*1e3,'--', 'Color',"#FFA500",'LineWidth',2); %#FFA500橙色
    l3 = xline(3* FWHM/2*1e3,'--', 'Color',"#FFD700",'LineWidth',2); %#FFD700金色

    % 图例和标题
    lgd = legend([l1,l2,l3],sprintf('FWHM = %.2f$\\mu$m',1000*FWHM), ...
        sprintf('2FWHM $\\eta_f = $%.2f \\%%',F(2* FWHM/2)*100), ...
        sprintf('3FWHM $\\eta_f = $%.2f \\%%',F(3* FWHM/2)*100));
    lgd.Interpreter = 'latex';
    lgd.FontSize = 16;
    lgd.FontName = "Times New Roman";
    lgd.Location = 'southeast';
    lgd.EdgeColor = '#D3D3D3';
    lgd.BackgroundAlpha = 0.6;
    xlabel('{\rho} ({\mu}m)',FontSize=16,FontName="Times New Roman");
    ylabel('Enclosed Energy Ratio',FontSize=16,FontName="Times New Roman");
    exportgraphics(f,fullfile(save_path,sprintf('enclosed_energy_ratio_f_%.1f.png', efl)),'Resolution',300)
    if show
        waitfor(f)
    else
        close(f);
    end

    %% 计算MTF
    [mtf_x, mtf_y] = calculate_MTF(I, G.step, 'zero_coef',1e5,'gpu_acceleration',gpu_acceleration);
    if show
        f = figure;
    else
        f = figure(Visible="off");
    end
    f.Position = [0 0 640 512];
    plot(G.axis_fft(numel(mtf_x)+1:end), mtf_x,'Color',"#4169E1",'LineWidth',2);  % #4169E1皇家蓝
    ax = gca;
    ax.FontSize = 16;
    ax.FontName = "Times New Roman";
    grid on
    hold on;
    plot(G.axis_fft(numel(mtf_y)+1:end), mtf_y,'--','Color',"#FFA500",'LineWidth',2); %#FFA500橙色
    % 图例和标题
    lgd = legend('{\it{x}}','{\it{y}}');
    lgd.FontSize = 16;
    lgd.FontName = "Times New Roman";
    lgd.Location = 'northeast';
    lgd.EdgeColor = '#D3D3D3';
    lgd.BackgroundAlpha = 0.6;
    xlabel('Frequency (lp/mm)',FontSize=16,FontName="Times New Roman");
    ylabel('Modulation Transfer Function',FontSize=16,FontName="Times New Roman");
    exportgraphics(f,fullfile(save_path,sprintf('MTF_f_%.1f.png', efl)),'Resolution',300)
    if show
        waitfor(f)
    else
        close(f);
    end
    
    %% 放大并插值光斑
    [x_interp,y_interp,I_interp] = extract_spot(I, G,center_x, center_y,12e-3,2);
    % 计算光强和中间值
    I = abs(I_interp).^2;
    mid_index_interp = round(length(y_interp) / 2);
    if show
        f = figure;
    else
        f = figure(Visible="off");
    end
    f.Position = [0 0 1280 448];
    % x-y截面
    subplot(1, 2, 1);   
    imagesc(x_interp*1e3, y_interp*1e3, I);
    ax=gca;
    ax.DataAspectRatio = [1,1,1];
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    ax.XLim = [-10 10];
    ax.YLim = [-10 10];
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Intensity';
    c.Label.FontSize = 16;
    c.Label.FontName = "Times New Roman";
    title('Interpolation Image',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{x}} ({\mu}m)',FontSize=16,FontName="Times New Roman");
    ylabel('{\it{y}} ({\mu}m)',FontSize=16,FontName="Times New Roman");
    % y 截面
    subplot(1, 2, 2);
    plot(x_interp*1e3, I(mid_index_interp, :),'Color',"#4169E1",'LineWidth',2); 
    ax=gca;
    ax.FontSize=16;
    ax.FontName = "Times New Roman";
    ax.XLim = [-10 10];
    title('{\it{y}}-axis Cross Section',FontSize=16,FontName="Times New Roman");
    xlabel('{\it{y}} ({\mu}m)',FontSize=16,FontName="Times New Roman");
    ylabel('Intensity',FontSize=16,FontName="Times New Roman");
    grid on;
    
    exportgraphics(f,fullfile(save_path,sprintf('image_interp_f_%.1f.png', efl)),'Resolution',300)

    if show
        waitfor(f)
    else
        close(f);
    end

    %% y-z截面计算
    if interval > 0 && sampling_point > 0
        disp('Calculating the light field in y-z cross section');
        dist_array = linspace(d_bfl - interval, d_bfl + interval, sampling_point);
        e_yz = zeros(length(G.axis), length(dist_array));
        progress_width = 50; % 进度条宽度     
        tic; % 记录起始时间
        output = 0;
        for i = 1:length(dist_array)
            iter_start = toc; % 记录单次循环开始时间
            p_bfl = PropOperator(G, S.wavelength_vacuum, dist_array(i), 'method', method, 'gpu_acceleration', gpu_acceleration);
            e_yz_d = p_bfl.prop(e_5);
            e_yz(:, i) = e_yz_d(mid_index_0, :);
            iter_time = toc - iter_start;
            remaining_time = iter_time*(length(dist_array)-i);
            percent_done = i / length(dist_array);
            num_hashes = round(percent_done * progress_width); % 进度条的数量
            fprintf(repmat('\b', 1, output)); % 删除百分比和时间字符
            output = fprintf('Progress: [%s%s] %3.0f%% | Iter Time: %.2fs | Time Remaining: %.2fs' ...
                ,repmat('=', 1, num_hashes),repmat(' ', 1, progress_width - num_hashes), ...
                percent_done * 100, iter_time, remaining_time);          
        end
    
        % 总耗时
        total_time = toc;
        fprintf('y-z cross section calculation is completed. Elapsed time: %.2f s\n', total_time);
    
        if show
            f = figure;
        else
            f = figure(Visible="off");
        end
        f.Position = [0 0 640 512];
        % x-y截面
        imagesc(dist_array, G.axis, log10(abs(e_yz).^2));
        ax=gca;
        ax.FontSize=16;
        ax.FontName = "Times New Roman";
        colormap(jet);
        c = colorbar;
        c.Label.String = 'log10(Intensity)';
        c.Label.FontSize = 16;
        c.Label.FontName = "Times New Roman";
        title('{\it{y}}-{\it{z}} Cross Section',FontSize=16,FontName="Times New Roman");
        xlabel('{\it{z}} (mm)',FontSize=16,FontName="Times New Roman");
        ylabel('{\it{y}} (mm)',FontSize=16,FontName="Times New Roman");
        exportgraphics(f,fullfile(save_path,sprintf('y-z_cross_section_f_%.1f.png', efl)),'Resolution',300)
        if show
            waitfor(f)
        else
            close(f);
        end

    end

       

end
    