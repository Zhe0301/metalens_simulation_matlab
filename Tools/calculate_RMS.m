function [centroid_x, centroid_y, rms_radius,total_energy] = calculate_RMS(I, G)
    % 计算质心和RMS半径
    % I: 光强矩阵
    % G: 坐标网路类
    I(I < max(I)/1e3)=0;
    total_energy = trapz(G.axis,trapz(G.axis,I,2)); %对光强进行坐标积分获得能量
    
    
    centroid_x = trapz(G.axis,trapz(G.axis,I .* G.d2_x,2)) / total_energy;% 计算质心位置
    centroid_y = trapz(G.axis,trapz(G.axis,I .* G.d2_y,2)) / total_energy;

    rms_radius = 2*sqrt(trapz(G.axis,trapz(G.axis,I .* ((G.d2_x - centroid_x).^2 + (G.d2_y - centroid_y).^2),2))/ total_energy);
end
