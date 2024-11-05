clear;
close all;
addpath '.\Tools';
% 配置保存路径
current_date = datetime('now', 'Format','yyyyMMdd');
name = 'actual_cylinder_1064';
save_path = sprintf('E:/Research/WavePropagation/metalens_simulation_matlab/Zoom_6x/%s_%s/', current_date, name);
if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf("文件夹 '%s' 已创建。\n", save_path);
else
    fprintf("文件夹 '%s' 已存在。\n", save_path);
end

% 建立仿真参数
period = 500e-6; % 单元结构周期，单位mm
grid_length = 3;
num_points = round(grid_length / period);
G = Grid(num_points, num_points * period);

% 结构参数
efl = [0.8, 1.6, 2.4, 3.2, 4, 4.8];
d_12 = [1.000000131063, 1.846692355066, 2.222548182242, 2.452026972933, 2.6206631695450, 2.776068441725];
d_23 = [3.442929821959, 2.341826841284, 1.798827702873, 1.444344104929, 1.186273959296, 0.9999998024761];
d_bfl = [2.307069517468, 2.561481226825, 2.728624214413, 2.853629959039, 2.943064336454, 2.973929452009];

% 波长和折射率
wavelength_vacuum = 1064 * 1e-6; % 单位mm
refractive_index = 1.5094756; % SiO2折射率

% 相位离散参数
discrete = 8;
unit_phase = [0, 0.771214, 1.611515, 2.329676, 3.129286, 3.941896, 4.700876, 5.486196];
unit_t = [0.957917, 0.885148, 0.799661, 0.756625, 0.738394, 0.713619, 0.756191, 0.81503];
points = linspace(0, 1, discrete + 1);
boundaries = arrayfun(@(i) [points(i), points(i + 1)], 1:discrete, 'UniformOutput', false);
boundaries = cell2mat(boundaries') * 2 * pi;

% 光源定义
s = Source(G, wavelength_vacuum, 1);
s.plane_wave(pi/2, pi/2);

% 镜片定义与仿真
L1 = Lens(G);
L1.binary2_d(0.72, [-5.848595615954e2, 3.20546768416e1, -1.422681034594e1, 2.984845643491], ...
    'm',1,'r_0',1, 'unit_phase', unit_phase, 'unit_t', unit_t, 'boundaries',boundaries);
% L1.plot_phase('show',true)

L2 = Lens(G);
L2.binary2_d(0.3, [3.157248960338e3, -2.062812665413e3, 1.207061495482e4, -6.173754564586e4], ...
    'm',1,'r_0',1, 'unit_phase', unit_phase, 'unit_t', unit_t, 'boundaries',boundaries);

L3 = Lens(G);
L3.binary2_d(1.2, [-1.845376080337e3, 6.935299381526e1, -7.928934283067, 1.711027879901], ...
    'm',1,'r_0',1, 'unit_phase', unit_phase, 'unit_t', unit_t, 'boundaries',boundaries);

% 执行三元变焦系统
d_lens = 0.75;
method = 'BL-AS';
for i = 1:length(efl)
    fprintf("EFFL = %.1f mm\n", efl(i));
    three_element_zoom_system_reverse(s, L1, L2, L3, G, d_lens, d_12(i), d_23(i), d_bfl(i), efl(i), ...
        'refractive_index',refractive_index, 'save_path',save_path,'sampling_point',50,'interval',1, ...
        'method', method, 'show',true, 'gpu_acceleration',true);
end
