function [skymask_fish, avail_ser_enu, pr_az_error] = BW_sn_Sky2skymask_upd(BW_sn_Sky, avail_ser_enu_oni)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk
% BW_SN_SKY2SKYMASK_UP 将二值天空图像转换为鱼眼天空掩码。
%   [skymask_fish, avail_ser_enu, pr_az_error] = BW_sn_Sky2skymask_up(BW_sn_Sky, avail_ser_enu_oni)
%   输入:
%     BW_sn_Sky         - 输入二值天空图像。
%     avail_ser_enu_oni - 可选的ENU序列（如果提供，重用）。
%   输出:
%     skymask_fish   - 鱼眼天空掩码（方位角和仰角）。
%     avail_ser_enu  - 处理后的ENU序列。
%     pr_az_error    - 方位角错误比例。

% 初始化和初步处理
[long, lat] = size(BW_sn_Sky);
I_line = sum(BW_sn_Sky, 1) ./ long;
Q_line = sum(BW_sn_Sky, 2);
Q_line_p = Q_line ./ lat;
hy_ele_id = find(diff(Q_line_p), 1, 'first');

% 调整上部掩码（清除上部框）
for idk = 1:lat
    fg_up = 0;
    fg_all_blk = 0;
    if all(BW_sn_Sky(1:hy_ele_id, idk) == 1)
        fg_up = 1;
    end
    if I_line(idk) < 0.99  % not all building
        fg_all_blk = 1;
    end
    if logical(fg_up .* fg_all_blk)
        BW_sn_Sky(1:hy_ele_id, idk) = 0.1;
    end
    
    lou_lab = find(BW_sn_Sky(:, idk) >= 1.5, 1);
    BW_sn_Sky(lou_lab:1:min([lou_lab + 5, long]), idk) = 1.5;
end

% 找到天空像素并转换为ENU坐标
[row_s, col_s, ~] = find(BW_sn_Sky >= 1);
backgroundInd = sub2ind(size(BW_sn_Sky), row_s, col_s);
avail_ser = [row_s, col_s, backgroundInd];
center_pic = [long, lat / 2];

if nargin <= 1
    avail_ser_enu_oni = [];
end

if ~isempty(avail_ser_enu_oni)
    avail_ser_enu = avail_ser_enu_oni;
    disp('ReUsing avail_ser_enu_oni')
else
    avail_ser_enu = [(avail_ser(:, 1:2) - center_pic), backgroundInd];
    az_sky = interp1([-center_pic(2), center_pic(2)], [0, 360], avail_ser_enu(:, 2));
    ele_sky = interp1([-center_pic(1), 0], [90, 0], avail_ser_enu(:, 1));
    avail_ser_enu(:, 5) = ele_sky;
    avail_ser_enu(:, 4) = roundn_anynum(az_sky, log10(0.2));
    avail_ser_enu = avail_ser_enu(~isnan(avail_ser_enu(:, 5)), :);
end

avail_ser_enu = unique(avail_ser_enu, "rows", "sorted");

% 唯一方位角和网格处理
uni_az = unique(avail_ser_enu(:, 4));
avail_ser_enu_N = round(avail_ser_enu);
[C, IA, IC] = unique(avail_ser_enu_N(:, 4:5), "rows", 'stable');
cc_grd_ORG = avail_ser_enu_N(IA, :);

% 计算占用率和参数
SKY_B = 5;
ocp_pre = sum(cc_grd_ORG(:, 5) >= (90 - SKY_B)) / SKY_B / 360;
ocp_pre(ocp_pre > 1) = 1;
dis = ocp_pre - [0, 0.159, 0.5, 1];
dis(dis > 0) = -inf;
[~, I] = max(dis);

a_low_con_1 = 0.8;
a_low_con = 2.2;

cc_grd=[cc_grd_ORG];
% 边界过滤循环（移除低置信度点）
exud_pre = 1;
while exud_pre > 0.005
    k = boundary(cc_grd(:, 4), cc_grd(:, 5), 1);
    bound_ele = [k, cc_grd(k, 4), cc_grd(k, 5)];
    
    az0_f = bound_ele(bound_ele(:, 2) == 0, :);
    [~, I0] = max(az0_f(:, 3));
    az360_f = bound_ele(bound_ele(:, 2) == 360, :);
    [~, I360] = max(az360_f(:, 3));
    
    bound_ele(bound_ele(:, 2) == 0, :) = [];
    bound_ele(bound_ele(:, 2) == 360, :) = [];
    bound_ele = bound_ele(bound_ele(:, 3) > 1, :);
    
    bound_ele = [bound_ele; az0_f(I0, :); az360_f(I360, :)];
    
    for id_p = 1:size(bound_ele, 1)
        exist_ele = cc_grd(cc_grd(:, 4) == bound_ele(id_p, 2), 5);
        exist_ele = exist_ele(exist_ele > 0);
        bound_ele(id_p, 4) = sum(exist_ele <= bound_ele(id_p, 3)) ./ bound_ele(id_p, 3);
        bound_ele(id_p, 5) = mean(exist_ele);
    end
    bound_ele = sortrows(bound_ele, 2);
    bound_ele(bound_ele(:, 3) < 0, :) = [];
    
    az_list = unique(bound_ele(:, 2));
    temp_bd = [];
    for ID = 1:length(az_list)
        temp_bd2 = [bound_ele];
        temp_bd2(temp_bd2(:, 2) ~= az_list(ID), 3) = nan;
        [M, I] = max(temp_bd2(:, 3));
        temp_bd = [temp_bd; bound_ele(I, :)];
    end
    bound_ele = temp_bd;
    
    low_con_1 = bound_ele(:, 4) < a_low_con_1;
    low_con = bound_ele(:, 3) > bound_ele(:, 5) .* a_low_con;
    low_con_2 = bound_ele(:, 3) > 90;
    low_con_a = [logical(low_con_1 + low_con + low_con_2)];
    
    cc_grd(bound_ele(low_con_a, 1), 5) = -1;
    exud_pre = sum(low_con_a) ./ length(bound_ele);
end

% 应用过滤结果
temp = cc_grd(IC, :);
avail_ser_enu(temp(:, 5) < 0, 4) = nan;
avail_ser_enu = avail_ser_enu(~isnan(avail_ser_enu(:, 4)), :);

% 仰角估计和错误计算
est_ele = [];
az_count = 0;
ele_biass = 2.5;
for id_az = 1:length(uni_az)
    cc_az_loc = find(avail_ser_enu(:, 4) == uni_az(id_az));
    cc_grd = avail_ser_enu(cc_az_loc, 4:5);
    Ls_dbsaca = 0;
    if size(cc_grd, 1) > 1
        cc_grd(:, 3:4) = [[99, 99]; diff(cc_grd(:, 1:2), 1)];
    else
        if isempty(cc_grd)
            cc_grd = [uni_az(id_az), 0.01];
        end
        cc_grd(:, 3:4) = [0, 0];
    end
    [M, ~] = min(cc_grd(:, 4));
    
    if M < -ele_biass
        az_count = az_count + 1;
        num_near_by = min([50, ceil(size(cc_grd, 1) * 0.1)]);
        idx_qt = 1;
        while length(unique(idx_qt)) <= 1
            idx_qt = dbscan(cc_grd(:, 1:2), ele_biass, ceil(num_near_by));
            num_near_by = num_near_by - 1;
            if num_near_by <= 2
                idx_qt = ones(size(cc_grd, 1), 1);
                break
            end
        end
        
        if any(idx_qt == -1)
            bad_ele = mean(avail_ser_enu(cc_az_loc(idx_qt == -1), 5));
            good_ele = sort((avail_ser_enu(cc_az_loc(idx_qt ~= -1), 5)), 'descend');
            good_ele_est = mean(good_ele(1:min([length(good_ele), 10])));
        else
            good_ele = sort((avail_ser_enu(cc_az_loc(idx_qt ~= -1), 5)), 'descend');
            good_ele_est = mean(good_ele(1:min([length(good_ele), 10])));
            bad_ele = max([-99, mean(avail_ser_enu(cc_az_loc(idx_qt == -1), 5))]);
        end
        
        if id_az == 1
            if any(idx_qt == -1)
                Ls_dbsaca = 1;
            end
        else
            Ls_dbsaca = abs(bad_ele - est_ele(id_az - 1, 2)) < abs(good_ele_est - est_ele(id_az - 1, 2));
        end
        
        if Ls_dbsaca == 1
            disp(['Ls_dbsaca- Processing--az' num2str(uni_az(id_az)), '-=P-' num2str(id_az / length(uni_az))])
            [avail_ser_enu_c, pick_idex] = az_range_ez(uni_az(id_az) + [-2, 2], avail_ser_enu(), 4);
            [lia, ~] = ismember(pick_idex, cc_az_loc, 'rows');
            avail_ser_enu_c = roundn_anynum(avail_ser_enu_c, log10(0.5));
            [C, IA, IC] = unique(avail_ser_enu_c(:, 4:5), "rows", 'stable');
            cc_grd = avail_ser_enu_c(IA, :);
            num_near_by = min([50, ceil(size(cc_grd, 1) * 0.1)]);
            idx_qt = 1;
            while length(unique(idx_qt)) <= 1
                idx_qt = dbscan(cc_grd(:, 4:5), ele_biass, ceil(num_near_by));
                num_near_by = num_near_by - 2;
                if num_near_by <= 2
                    idx_qt = ones(size(cc_grd, 1), 1);
                    break
                end
            end
            idx_qt = idx_qt(IC);
            idx_qt = idx_qt(lia);
        end
        
        avail_ser_enu(cc_az_loc(idx_qt == -1), 4) = nan;
        good_ele_est = max(avail_ser_enu(cc_az_loc(idx_qt ~= -1), 5));
        est_ele(id_az, 1:2) = [uni_az(id_az), good_ele_est];
    else
        est_ele(id_az, 1:2) = [uni_az(id_az), max(cc_grd(:, 2))];
    end
end

pr_az_error = az_count ./ length(uni_az);
avail_ser_enu = avail_ser_enu(~isnan(avail_ser_enu(:, 4)), :);
avail_ser_enu(:, 4) = mod(avail_ser_enu(:, 4) + 0, 360);

% 拟合天空掩码
cluster_max_p = max_sat_each    (avail_ser_enu(:, 5), avail_ser_enu(:, 4));
gap = 1;
cluster_max_p(:, 5) = roundn_anynum(cluster_max_p(:, 1), log10(gap));
skymask_fish_in = nan(360 ./ gap + 1, 2);
skymask_fish_in(:, 1) = 0:gap:360;

for idx = 1:length(skymask_fish_in)
    mean_bu = mean(cluster_max_p(cluster_max_p(:, 5) == skymask_fish_in(idx, 1), 2));
    if ~isempty(mean_bu)
        skymask_fish_in(idx, 2) = mean_bu;
    end
end

skymask_fish_in(isnan(skymask_fish_in(:, 2)), :) = [];

sk2ft=[[skymask_fish_in(:, 1)-360;skymask_fish_in(:, 1);skymask_fish_in(:, 1)+360], ...
    [skymask_fish_in(:, 2);skymask_fish_in(:, 2);skymask_fish_in(:, 2)]];

sk2ft=unique(sk2ft,"rows");

f0 = fit(sk2ft(:,1), sk2ft(:,2),'linearinterp');

BU_EST = feval(f0, 0:1:360);
skymask_fish = [];
skymask_fish(:, 1) = 0:1:360;
skymask_fish(:, 2) = BU_EST;

if skymask_fish(1, 2) ~= skymask_fish(end, 2)
    mean_0_360 = mean(mean(skymask_fish([1, end], 2)));
    skymask_fish([1, end], 2) = mean_0_360;
end

if ocp_pre > 0.1
    pr_az_error = -pr_az_error;
    skymask_fish([1, ], 2) = nan;
end

end
