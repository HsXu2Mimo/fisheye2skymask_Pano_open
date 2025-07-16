function [com_2W_w_up, restart_fg] = gap_in_seg_rep_upd(com_2W_w, Half_J_OUT)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk
% GAP_IN_SEG_REP 修复图像分割中的间隙。
%   [com_2W_w_up, restart_fg] = gap_in_seg_rep(com_2W_w, Half_J_OUT)
%   输入:
%     com_2W_w    - 输入分割矩阵。
%     Half_J_OUT  - 半图像输出。
%   输出:
%     com_2W_w_up - 更新后的分割矩阵。
%     restart_fg  - 重启标志（0或1）。

% 初始化
restart_fg = 0;
[long, lat] = size(com_2W_w);
Cap_pre = [];
Cap_pre_y_Max = [];

% 计算每个列的容量和最大Y值
for i1dx = 1:lat
    lin_0 = find(com_2W_w(:, i1dx) >= 1);
    lin_2 = find(com_2W_w(:, i1dx) == 2.5);
    lin_3 = max([nan; find(com_2W_w(:, i1dx) == 1.5)]);
    
    Cap_pre_y_Max(i1dx, 1) = lin_3;
    Cap_pre(i1dx, 1) = (length(lin_0) - length(lin_2)) ./ long;
end

% 拟合和过滤容量值
Cap_pre_fb = sort(Cap_pre);
[mu, s, muci, sci] = normfit(Cap_pre_fb);
Cap_pre(Cap_pre < max([0.1, mu + s .* 1.5])) = nan;

% 使用kmedoids聚类
try
    Cap_pre = reshape(Cap_pre, [], 1);
    [iq_L, C, SUMD] = kmedoids([(1:length(Cap_pre))', reshape(Cap_pre, [], 1)], 2);
    iq = [(1:length(Cap_pre))', reshape(Cap_pre, [], 1), iq_L];
catch
    iq = 3 .* ones(size(Cap_pre, 1), 4);
end

% 初始化更新矩阵
com_2W_w_up = sign(com_2W_w);

% 循环处理每个聚类
for id_iq = 1:2
    % 定义更新范围
    gap_update_1 = max([1, min(iq(iq(:, 3) == id_iq, 1)) - 70]):1:...
        min([max(iq(iq(:, 3) == id_iq, 1)) + 70, lat]);
    
    gap_update_2 = max([1, min(iq(iq(:, 3) == id_iq, 1) + 10)]):1:...
        min([max(iq(iq(:, 3) == id_iq, 1) - 10), lat]);
    
    Y = long;
    
    % 处理第一个范围
    com_2W_w_bOTH = zeros(size(com_2W_w));
    BW_hf_ec = Half_J_OUT(1:Y, gap_update_1, :);
    com_2W_w_p = com_2W_w(1:Y, gap_update_1);
    Lout_00t{1} = Bn_lout_SMALL(BW_hf_ec, com_2W_w_p);
    
    % 处理第二个范围（如果足够大）
    if length(gap_update_2) > 200
        BW_hf_ec = Half_J_OUT(:, gap_update_2, :);
        com_2W_w_p = com_2W_w(1:Y, gap_update_2);
        Lout_00t{2} = Bn_lout_SMALL(BW_hf_ec, com_2W_w_p);
    else
        Lout_00t{2} = Lout_00t{1};
        gap_update_2 = gap_update_1;
    end
    [L, Centers] = imsegkmeans( single( rgb2gray(BW_hf_ec)),2, 'NumAttempts',6);

    [bu_Col,cand_color,top_color]=Sky_bu_col(L);
    Lout_00t{3}=L~=top_color;

    % 更新com_2W_w_bOTH
    com_2W_w_bOTH(1:Y, gap_update_1) = com_2W_w_bOTH(1:Y, gap_update_1) + Lout_00t{1};
    com_2W_w_bOTH(1:Y, gap_update_2) = com_2W_w_bOTH(1:Y, gap_update_2) + ...
        sign(Lout_00t{3} + Lout_00t{2});
    
    % 二值化和前景/背景分割
    Bw_com = double(com_2W_w_bOTH(:, gap_update_2) >= 1);
    [BW_sn_Sky_1] = segem_fore_back(BW_hf_ec, Bw_com);
    [BW_sn_Sky_2] = segem_fore_back(BW_hf_ec, Bw_com);
    BW_sn_Sky = BW_sn_Sky_1 + BW_sn_Sky_2;
    
    % 更新输出
    com_2W_w_up(:, gap_update_2) = double(BW_sn_Sky) + (Bw_com);
end

end

% 子函数：小区域分割
function Lout_00t6 = Bn_lout_SMALL(BW_hf_ec, com_2W_w_p)
% BN_LOUT_SMALL 小区域的颜色分割和处理。
%   Lout_00t6 = Bn_lout_SMALL(BW_hf_ec, com_2W_w_p)

pre_sk = inf(6, 3);
for color_num = 3:6
    [Lout_001, ~, ~] = rgb_BW2imseg((BW_hf_ec), color_num);
    Lout_00 = Lout_001;
    [~, ~, top_color] = Sky_bu_col(Lout_001);
    
    [L_m0, L_n0] = size(Lout_00);
    prec_sky(color_num) = sum(Lout_00 == top_color, 'all') ./ L_m0 ./ L_n0;
    
    CC_bu_PRE = [];
    CC_bu_PRE_s = [];
    for cc_cl = 1:color_num
        CC_bu_PRE(cc_cl) = sum(sign(com_2W_w_p) .* (Lout_00 == cc_cl), 'all') ./ sum((Lout_00 == cc_cl), 'all');
        CC_bu_PRE_s(cc_cl) = sum(sign(com_2W_w_p) .* (Lout_00 == cc_cl), 'all') ./ numel(com_2W_w_p);
    end
    
    sky_bas = max([0.25, min(CC_bu_PRE)]);
    
    Lout_00t = zeros(size(Lout_00));
    Lout_00t(ismember(Lout_00, find(CC_bu_PRE() > sky_bas))) = 1;
    Lout_00t_C{color_num} = Lout_00t;
    
    pre_sk(color_num, 1:2) = [sum(CC_bu_PRE() > sky_bas) ./ color_num, ...
        sum(CC_bu_PRE_s(CC_bu_PRE > sky_bas))];
    pre_sk(color_num, 1 + 2) = pre_sk(color_num, 1) .* pre_sk(color_num, 2);
end

[~, I] = min(pre_sk(:, 3));
Lout_00t6 = Lout_00t_C{I};

end

% 旧版本代码（作为参考保留，注释掉）
% function old_d
% % OLD_D 旧版本的间隙修复逻辑。
% for id_iq = 1:2
%     gap_update = max([1, min(iq(iq(:,4) == id_iq,1)) - 70]):1:...
%         min([max(iq(iq(:,4) == id_iq,1)) + 70, lat]);
%     
%     BW_hf_ec = Half_J_OUT(:, gap_update, :);
%     com_2W_w_p = com_2W_w(:, gap_update);
%     if length(gap_update) > 200
%         restart_fg = 1;
%         prec_sky = [];
%         
%         for color_num = 3:6
%             [Lout_001, ~] = rgb_BW2imseg((BW_hf_ec), color_num);
%             [bu_Col, cand_color, top_color] = Sky_bu_col(Lout_001);
%             Lout_00 = Lout_001;
%             [L_m0, L_n0] = size(Lout_00);
%             prec_sky(color_num) = sum(Lout_00 == top_color, 'all') ./ L_m0 ./ L_n0;
%             if abs(prec_sky(color_num) - prec_sky(color_num - 1)) < 0.1
%                 CC_bu_PRE = [];
%                 for cc_cl = 1:color_num
%                     CC_bu_PRE(cc_cl) = sum(sign(com_2W_w_p) .* (Lout_00 == cc_cl), 'all') ./ sum((Lout_00 == cc_cl), 'all');
%                 end
%                 sky_bas = max([0.25, min(CC_bu_PRE)]);
%                 for cc_cl = 1:color_num
%                     if CC_bu_PRE(cc_cl) > sky_bas
%                         Lout_00(Lout_00 == cc_cl) = -top_color;
%                     else
%                         Lout_00(Lout_00 == cc_cl) = top_color;
%                     end
%                 end
%                 break;
%             end
%         end
%         
%         [BW_sn_Sky_1] = segem_fore_back(BW_hf_ec, Lout_00 ~= top_color);
%         [BW_sn_Sky_2] = segem_fore_back(BW_hf_ec, Lout_00 ~= top_color);
%         BW_sn_Sky = BW_sn_Sky_1 + BW_sn_Sky_2;
%         
%         com_2W_w_up(:, gap_update) = double(BW_sn_Sky) + (Lout_00 ~= top_color);
%     end
% end
% restart_fg = 0;
% end
