function [skymask_fish, skymask_image, Corr_mark] = fisheye2skymask_Pano_up_fix2hwSP...
    (pic_name, heading_rotation, corrtion_1, plot_flag)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

% 初始化随机数生成器和变量
rng('default')
Corr_mark = 1;
restart_fg = 0;
skymask_image = [];
disp(['Processing the Panorama2skyplot imAge // '])
disp([pic_name])

% 读取图像并进行初始转换
[pathstr, name, ext] = fileparts(pic_name);
I = imread(pic_name);
clor_resloution = 30 * 1024;
[CW, map] = rgb2ind(I, clor_resloution);
[J_OUT, J_seg] = img_I2BW_upd(I);
[w_i, h_i, ~] = size(J_OUT);

% 准备分块参数
mid_pic = linspace(0, h_i, 3);
med_pic = cumsum(linspace(0, h_i, 3)) ./ 2;
mid_pic_List = [[0, med_pic(end)]; [med_pic(2), mid_pic(end)]];
mid_pic_List(:, 1) = mid_pic_List(:, 1) + 1;

% 初始化组合掩码矩阵
com_2W_w = zeros(w_i/2, h_i);
com_2W_b = zeros(w_i/2, h_i);
com_2W_L = zeros(w_i/2, h_i);
BW_hf_C = [];
BW_hf_C_SEG = [];

% 预处理上半部分图像
Half_J_OUT = J_OUT(1:end/2, :, :);
if nargin > 3
    [Lout, ~] = rgb_BW2imseg_upd(Half_J_OUT, 2);
else
    [Lout, ~] = rgb_BW2imseg_upd(Half_J_OUT);
end
[bu_Col, cand_color, top_color] = Sky_bu_col(Lout);

stg_color = cand_color(cand_color(:, 4) <= 0.05, 1);
if ~isempty(stg_color)
    disp(['Too Obvious color with ' num2str(cand_color(cand_color(:, 4) <= 0.05, 4))])
    [row_s, col_s, ~] = find(Lout ~= stg_color);
    backgroundInd = sub2ind(size(Lout), row_s, col_s);
    [Alab, map_obs] = rgb2ind(Half_J_OUT, clor_resloution);
    perm_n = cumprod(size(Alab));
    perm_n = 1:perm_n(2);
    perm_n(backgroundInd) = [];
    Alab(perm_n) = nan;
    Half_J_OUT = ind2rgb(Alab, map_obs);

    [Alab, map_obs] = rgb2ind(J_seg(1:end/2, :, :), clor_resloution);
    Alab(perm_n) = nan;
    Half_J_seg = ind2rgb(Alab, map_obs);
else
    Half_J_seg = J_seg(1:end/2, :, :);
end



% 分块处理上半部分图像
for id_block = 1:size(mid_pic_List, 1)
    BW_hf_ec = Half_J_OUT(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2), :);
    BW_hf_ec_SEG = Half_J_seg(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2), :);

    [Lout_0, ~] = rgb_BW2imseg_upd(BW_hf_ec, 3);
    [~, cand_color, ~] = Sky_bu_col(Lout_0);
    ex_color = cand_color(cand_color(:, 4) <= 0.05, 1);
    BW_hf_ec = img_Lout_clear(BW_hf_ec, Lout_0, ex_color);%%$ clean the unnec-col
    [Lout, ~, ~] = rgb_BW2imseg_upd(BW_hf_ec, 3);
    % figure
    % imagesc(BW_hf_ec)
    [bu_Col, cand_color, top_color] = Sky_bu_col(Lout);

    [BW_sn_Sky_1] = segem_fore_back_upd(BW_hf_ec_SEG, Lout ~= top_color);
    [BW_sn_Sky_2] = segem_fore_back_upd(BW_hf_ec_SEG, Lout ~= top_color);

    BW_sn_Sky = BW_sn_Sky_1 + BW_sn_Sky_2;
    com_2W_b(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2)) = ...
        com_2W_b(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2)) + BW_sn_Sky;
    BW_hf_C(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2), :) = BW_hf_ec;
    % BW_hf_C_SEG(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2), :) = BW_hf_ec_SEG;
    com_2W_L(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2)) = ...
        com_2W_L(:, mid_pic_List(id_block, 1):mid_pic_List(id_block, 2)) + double(Lout ~= top_color);

end

% 生成初始组合掩码并检查差距
com_2W_w0 = sign(com_2W_b) + 1.5 .* sign(com_2W_L);
% [com_2W_w, restart_fg] = gap_in_seg_rep_upd(com_2W_w0, Half_J_OUT);
bw = activecontour(Half_J_OUT, com_2W_w0==0, 300, 'Chan-Vese',1);
% 生成天空掩码
[skymask_fish, avail_ser_enu, pr_az_error] = BW_sn_Sky2skymask_upd(bw==0);


% 可选重启绘图（如果 plot_flag 存在）
if nargin > 3
    figure
    subplot(4, 2, 1:2)
    imagesc(I(1:end/2, :, :))
    set(gca, 'Colormap', map);
    title(['ST1\', name])
    xticks([0:200:h_i-100, h_i])

    subplot(4, 2, 3:4)
    imagesc(sign(com_2W_L))
    xticks([0:300:h_i-100, h_i])

    % subplot(4, 2, 4)
    % imagesc(sign(bw==0))
    % xticks([0:300:h_i-100, h_i])

    subplot(4, 2, 5:6)
    imagesc(bw==0)
    xticks([0:200:h_i-100, h_i])

    subplot(4, 2, 7:8)
    [Alab, map_sb] = rgb2ind(Half_J_OUT, clor_resloution);
    set(gca, 'Colormap', map_sb);
    coco = double(Alab(avail_ser_enu(:, 3)));
    scatter(avail_ser_enu(:, 4), avail_ser_enu(:, 5), 4, map_sb(double(coco) + 1, :), 'filled')
    % plot_az_label;
    hold on
    xlim([0,360])
    plot(skymask_fish(:, 1), skymask_fish(:, 2), 'r-', 'LineWidth', 2)
    set(gcf, 'position', [76 088 800 900]);
end
% 生成最终输出（应用旋转和修正）
if nargin > 1
    skymask_image = skymask_transfer(skymask_fish, heading_rotation, 0);
end
if nargin > 2
    skymask_image = skymask_transfer(skymask_fish, 0, corrtion_1);
end