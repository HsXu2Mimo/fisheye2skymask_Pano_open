function [BW_sn_Sky, I_BW, foresub, backgroundsub, OK_TAG, foresub_bk] = ....
    segem_fore_back_upd(BW_hf_ec, Lout, plot_flg)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

% SEGEM_FORE_BACK 使用Lazy Snapping进行前景/背景分割。
%   [BW_sn_Sky, I_BW, foresub, backgroundsub, OK_TAG, foresub_bk] = segem_fore_back(BW_hf_ec, Lout, plot_flg)
%   输入:
%     BW_hf_ec  - 输入图像。
%     Lout      - 分割标签矩阵（可选）。
%     plot_flg  - 绘图标志（如果提供，显示结果）。
%   输出:
%     BW_sn_Sky      - 分割后的二值图像（天空/背景）。
%     I_BW           - 处理后的图像。
%     foresub         - 前景子采样点。
%     backgroundsub  - 背景子采样点。
%     OK_TAG         - 成功标签（1为成功，或调整后的值）。
%     foresub_bk     - 备份前景子采样点。

% 初始化
OK_TAG = 1;

% 图像预处理：如果未提供Lout，进行灰度转换和S曲线调整
if nargin < 2
    I_BW = rgb2gray(BW_hf_ec);
    [h, w] = size(I_BW);
    gray = I_BW;
    NewImage1 = zeros(h, w);
    a = 80 / 256;
    b = 180 / 256;
    c = 30 / 256;
    d = 220 / 256;
    for x = 1:w
        for y = 1:h
            if gray(y, x) < a
                NewImage1(y, x) = gray(y, x) * c / a;
            elseif gray(y, x) < b
                NewImage1(y, x) = (gray(y, x) - a) * (d - c) / (b - a) + c;
            else
                NewImage1(y, x) = (gray(y, x) - b) * (255 - d) / (255 - b) + d;
            end
        end
    end
    I_BW = NewImage1;
else
    I_BW = (BW_hf_ec);
end

% 获取图像尺寸
h = size(I_BW, 1);
w = size(I_BW, 2);
try
    [L_m, L_n] = size(Lout);
catch
    [L_m, L_n] = size(I_BW);
end

% 初始化前景和背景子采样
foresub = [];
backgroundsub = [];
foresub_bk = [nan, nan];

% 计算主颜色和非主颜色比例
numColors = length(unique(Lout));
cand_color = [];
cand_color(:, 1) = [unique(Lout)];
for co_id = 1:numColors
    [i, j] = find(Lout == cand_color(co_id, 1));
    cand_color(co_id, 2:3) = mean([i, j], 1);
end
[~, I] = max(cand_color(:, 2));
master_color = cand_color(I, 1);

counter_bu = (sum(Lout ~= master_color, 'all') ./ L_m ./ L_n);
disp(['?!!Building Inside ->-' , '->Color-Mode->' num2str(counter_bu)]);

% 循环遍历列，生成前景和背景标记
Ac = 0;
for i1dx = 1:L_n
    lin_0 = find(Lout(:, i1dx) == master_color);
    lin_02 = find(Lout(:, i1dx) ~= master_color);
    lin_0 = rmoutliers(lin_0, 'grubbs');
    
    if ~isempty(lin_0)
        lin_1 = lin_0(1);
        lin_1 = lin_1 - 1;
    else
        lin_1 = nan;
    end
    
    if lin_1 <= 0
        lin_1 = 1;
    end
    
    fr_l = 1:10:lin_1;
    
    if abs(length(lin_0) ./ L_m - (L_m - lin_1) ./ L_m) > 0.15
        ran_k = randperm(length(lin_02), ceil(min([length(lin_02) ./ 5, 20])));
        Ac = Ac + 1;
        foresub_bk = [foresub_bk; [reshape(fr_l, [], 1), ones(length(fr_l), 1) .* i1dx]];
    end
    
    foresub = [foresub; [reshape(fr_l, [], 1), ones(length(fr_l), 1) .* i1dx]];
    
    if ~isempty(lin_0)
        bk_l = max([lin_1, fix(h - h ./ 10)]):10:L_m;
    else
        bk_l = nan;
    end
    
    backgroundsub = [backgroundsub; [reshape(bk_l, [], 1), ones(length(bk_l), 1) .* i1dx]];
end

% 处理NaN和随机采样
format longG

foresub = foresub(~isnan(foresub(:, 1)), :);
fore_rand = randperm(size(foresub, 1), ceil(min([size(foresub, 1) ./ 5, 0.3e4])));

backgroundsub = backgroundsub(~isnan(backgroundsub(:, 1)), :);
Bg_rand = randperm(size(backgroundsub, 1), ceil(min([size(backgroundsub, 1) ./ 2, 0.3e4])));

% 生成前景和背景索引
fore_sub_n_nan = [[foresub(fore_rand, 1); foresub_bk(:, 1)], [foresub(fore_rand, 2); foresub_bk(:, 2)]];
fore_sub_n_nan(any(isnan(fore_sub_n_nan), 2), :) = [];

foregroundInd = sub2ind([h, w], fore_sub_n_nan(:, 1), fore_sub_n_nan(:, 2));
foregroundInd = unique(foregroundInd);
foregroundInd = foregroundInd(~isnan(foregroundInd(:, 1)), :);

backgroundInd = sub2ind([h, w], backgroundsub(Bg_rand, 1), backgroundsub(Bg_rand, 2));

% 计算superpixels
C_S = cumprod(size(I_BW));
L = superpixels(I_BW, ceil(C_S(2) / 5));

% 执行Lazy Snapping
% tic
disp(['Emty pre-->' num2str(Ac ./ L_n)])
BW_sn_Sky = lazysnapping(I_BW, L, foregroundInd, backgroundInd, ...
    'Connectivity', 8, 'EdgeWeightScaleFactor', 1000);
% toc

% 额外灰度版本处理
L = superpixels(rgb2gray(BW_hf_ec), ceil(C_S(2) / 5));
BW_sn_Sky_grey = lazysnapping(rgb2gray(BW_hf_ec), L, foregroundInd, backgroundInd, ...
    'Connectivity', 8, 'EdgeWeightScaleFactor', 1000);
BW_sn_Sky = BW_sn_Sky + BW_sn_Sky_grey;

% 后处理
BW_sn_Sky = BW_sn_Sky == 0;
BW_sn_Sky = double(BW_sn_Sky);

% 可选绘图
if nargin > 2
    figure;
    imagesc(BW_sn_Sky);
    set(gcf, 'position', [100, 100, 500, 200]);
end

% 如果counter_bu过高，调整输出
if counter_bu > 0.9
    disp(['?!!Building Inside -Too MUCH ending'])
    OK_TAG = BW_sn_Sky;
    BW_sn_Sky = master_color .* ones(h, w);
end

end
