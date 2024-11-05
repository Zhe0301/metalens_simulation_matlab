function [center_x, center_y, r_x_p, r_x_n, r_y_p, r_y_n] = calculate_FWHM(I, G,varargin)       
% 计算光强最强位置和半高全宽FWHM
% 参数:
% I: 光强矩阵
% G: 坐标网路类
% 可选参数:
% centrosymmetry: 是否中心对称
% 返回值:
% center_x, center_y: 最强光强位置
% r_x_p, r_x_n, r_y_p, r_y_n: x,y正负方向半高全宽半径
    p = inputParser;            
    addParameter(p,'centrosymmetry',true);      
    parse(p,varargin{:}); 
    centrosymmetry = p.Results.centrosymmetry;
    [max_val, max_idx] = max(I(:)); 
    [max_row, max_col] = ind2sub(size(I), max_idx);
    
    max_all_idx = find(abs(I-max_val)<=max_val*1e-5); % 
    if length(max_all_idx)>1
        center_x = sum(G.d2_x(max_all_idx))/length(max_all_idx);
        center_y = sum(G.d2_y(max_all_idx))/length(max_all_idx);
        
    else
        center_x = G.d2_x(max_row,max_col);
        center_y = G.d2_y(max_row,max_col);
    end
    max_val = max_val/2; % 光强的一半
    if centrosymmetry && center_x == center_y    
        I_x_p = I(max_row,max_col:end);
        x_p = G.d2_x(max_row,max_col:end);
        I_x_n = flip(I(max_row,1:max_col));
        x_n = flip(G.d2_x(max_row,1:max_col));
        I_y_p = I(max_row:end,max_col);
        y_p = G.d2_y(max_row:end,max_col);
        I_y_n = flip(I(1:max_row,max_col));
        y_n = flip(G.d2_y(1:max_row,max_col));
    
        position = find(I_x_p<max_val/2,1); % 小于最大光强一半的第一个值，x正方向
        r_x_p=interp1(I_x_p(1:position),x_p(1:position),max_val)-center_x;
        position = find(I_x_n<max_val/2,1); % x负方向
        r_x_n=-interp1(I_x_n(1:position),x_n(1:position),max_val)+center_x;
        position = find(I_x_p<max_val/2,1); % y正方向
        r_y_p=interp1(I_y_p(1:position),y_p(1:position),max_val)-center_y;
        position = find(I_x_p<max_val/2,1); % y负方向
        r_y_n=-interp1(I_y_n(1:position),y_n(1:position),max_val)+center_y;
    else
        I_x_p = I(max_row,max_col:end);
        x_p = G.d2_x(max_row,max_col:end);
        position = find(I_x_p<max_val/2,1); % 小于最大光强一半的第一个值
        r_x_p=interp1(I_x_p(1:position),x_p(1:position),max_val)-center_x;
        r_x_n=r_x_p; % 中心对称全都相等
        r_y_p=r_x_p;
        r_y_n=r_x_p;
    end
    

  
end
