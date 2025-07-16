function [Lout, bu_Col, Lout_2C] = rgb_BW2imseg_upd(BW_hf, numColors, cloplot_flg)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk
% RGB_BW2IMSEG 使用superpixels和kmeans进行图像分割。
%   [Lout, bu_Col, Lout_2C] = rgb_BW2imseg(BW_hf, numColors, cloplot_flg)
%   输入:
%     BW_hf       - 输入RGB图像。
%     numColors   - 聚类颜色数量（默认2）。
%     cloplot_flg - 绘图标志（如果提供，显示图像和标签）。
%   输出:
%     Lout     - 分割标签矩阵。
%     bu_Col   - 建筑物颜色ID。
%     Lout_2C  - 二值化标签矩阵。

% 初始化输出
bu_Col = [];
Lout_2C = [];

% 处理默认参数
if nargin == 1
    numColors = 2;
end

% 设置颜色分辨率并转换到LAB空间
clor_resloution = 30 * 1024;
Alab = rgb2lab(BW_hf);

% 获取图像尺寸
[m1, m2, ~] = size(Alab);

% Superpixels计算和均值颜色处理（循环处理superpixel数量）
for superpixel_num = ceil(m1 .* m2)
    % 计算superpixels
    [L, N] = superpixels(Alab, superpixel_num, 'isInputLab', true);
    
    % 获取每个superpixel的像素索引
    pixelIdxList = label2idx(L);
    [m, n] = size(L);
    
    % 计算每个superpixel的均值颜色
    meanColor = zeros(m, n, 3, "single");
    for i = 1:N
        meanColor(pixelIdxList{i}) = mean(Alab(pixelIdxList{i}));
        meanColor(pixelIdxList{i} + m * n) = mean(Alab(pixelIdxList{i} + m * n));
        meanColor(pixelIdxList{i} + 2 * m * n) = mean(Alab(pixelIdxList{i} + 2 * m * n));
    end
    
    % 使用kmeans进行分割
    [Lout, Centers] = imsegkmeans(meanColor, numColors, 'numAttempts', 6);
    
    % 计算主颜色（顶部50行的平均标签）
    [L_m, L_n] = size(Lout);
    master_color = round(mean(Lout(1:50, :), 'all'));
    
    % 计算非主颜色的比例
    counter_bu = (sum(Lout ~= master_color, 'all') ./ L_m ./ L_n);
end

% 可选绘图（如果cloplot_flg提供）
if nargin > 2
    figure;
    subplot(2, 1, 1)
    imagesc(Lout);
    subplot(2, 1, 2)
    imagesc(BW_hf);
end

% 处理颜色候选并确定建筑物颜色（如果numColors >= 2）
if numColors >= 2
    cand_color = [];
    cand_color(:, 1) = [1:numColors];
    for co_id = 1:numColors
        [i, j] = find(Lout == co_id);
        if isempty(i)
            i = nan;
            j = nan;
        end
        cand_color(co_id, 2:3) = mean([i, j], 1);
    end
    [~, I] = max(cand_color(:, 2));
    bu_Col = cand_color(I, 1);
    Lout_2C = Lout;
    Lout_2C(Lout_2C ~= bu_Col) = 0;
end

% 根据counter_bu调整输出并显示信息
if counter_bu > 0.2
    disp(['?!!Building Inside ->superpixel_num-' num2str(superpixel_num) ...
        , '->Color-Mode->' num2str(counter_bu)])
else
    disp(['?!!Building Inside ->superpixel_num-' num2str(superpixel_num) ...
        , '->Color-Mode->' num2str(counter_bu)])
    if nargin > 2
        figure;
        subplot(2, 1, 1)
        imagesc(Lout);
        subplot(2, 1, 2)
        imagesc(BW_hf);
    end
    Lout = Lout_2C;
end

% 转换为double类型
Lout_2C = double(Lout_2C);
Lout = double(Lout);

end
