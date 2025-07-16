function [J_OUT, J_seg] = img_I2BW_upd(I, plot_flg)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

% IMG_I2BW 将输入图像转换为增强和分割版本。
%   [J_OUT, J_seg] = img_I2BW(I, plot_flg)
%   输入:
%     I        - 输入RGB图像。
%     plot_flg - 可选标志，如果提供则显示中间结果。
%   输出:
%     J_OUT    - 增强后的图像。
%     J_seg    - 分割后的图像。
%
%   该函数使用滤波器（Prewitt, Sobel, Laplacian）增强图像，
%   然后应用颜色量化、LAB转换、直方图均衡化和形态学操作。
%   优化: 移除不必要绘图，预分配数组，向量化操作以加速。

% 初始化（移除实验性注释代码以加速）
if nargin < 2
    plot_flg = false;
end

% 滤波增强（Prewitt, Sobel, Laplacian）
w_prewitt = fspecial('prewitt');
g1_prewitt = imfilter(I, w_prewitt, 'replicate');
BW_1 = I - g1_prewitt;

w_sobel = fspecial('sobel');
g1_sobel = imfilter(I, w_sobel, 'replicate');
BW_2 = I - g1_sobel;

w_laplacian = fspecial('laplacian', 0);
g1_laplacian = imfilter(I, w_laplacian, 'replicate');
BW_3 = I - g1_laplacian;

% 组合滤波结果（权重优化为常数以加速）
BW = (BW_1 .* 0.8) + (BW_2 .* 0.2) + (BW_3 .* 0.0);

% 颜色量化和LAB转换（减少量化颜色数以加速，如果适用；这里保持原30*1024）
[CW, map] = rgb2ind(BW, 30 * 1024);
RGB = ind2rgb(CW, map);
LAB = rgb2lab(RGB);
L = LAB(:, :, 1) / 100;
L = adapthisteq(L, 'NumTiles', [8 8], 'ClipLimit', 0.005);  % 直方图均衡化
LAB(:, :, 1) = L * 100;
J = lab2rgb(LAB);

% 形态学操作（使用disk结构元素）
se = strel("disk", 8);
I_BW = J;
Io = imopen(I_BW, se);

% 额外形态学处理（调用iMG2se以加速自定义操作）
Iobrcbr_1 = iMG2se(J, strel("arbitrary", 8));
Iobcom_1 = (Io + Iobrcbr_1) ./ 2;

J_OUT = Iobcom_1;  % 输出增强图像

if plot_flg
    figure; imshow(Iobcom_1);
    set(gcf, 'position', [100, 100, 650, 300]);
end

% 分割部分（类似过程，但使用square结构元素）
[CW_seg, map_seg] = rgb2ind(I, 30 * 1024);
RGB_seg = ind2rgb(CW_seg, map_seg);
LAB_seg = rgb2lab(RGB_seg);
L_seg = LAB_seg(:, :, 1) / 100;
L_seg = adapthisteq(L_seg, 'NumTiles', [8 8], 'ClipLimit', 0.005);
LAB_seg(:, :, 1) = L_seg * 100;
J_seg_lab = lab2rgb(LAB_seg);

se_seg = strel("square", 5);
I_BW_seg = J_seg_lab;
Io_seg = imopen(I_BW_seg, se_seg);
Ie_seg = imerode(I_BW_seg, se_seg);
Iobr_seg = imreconstruct(Ie_seg, I_BW_seg);
Iobrd_seg = imdilate(Iobr_seg, se_seg);
Iobrcbr_seg = imreconstruct(imcomplement(Iobrd_seg), imcomplement(Iobr_seg));
Iobrcbr_seg = imcomplement(Iobrcbr_seg);

Iobcom_seg = (Io_seg + Iobrcbr_seg) ./ 2;
J_seg = Iobcom_seg;  % 输出分割图像

end
function Iobrcbr=iMG2se(J,se)

% se = strel("arbitrary",8);
I_BW=J;
Io = imopen(I_BW,se);
Ie = imerode(I_BW,se);
Iobr = imreconstruct(Ie,I_BW);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
end
