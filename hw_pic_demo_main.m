
% Copyright (C) 2020-2025 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

clc
clear all
close all
Pic_name='demo.jpg';


[skymask_fish,~,Corr_mark]=....
    fisheye2skymask_Pano_up_fix2hwSP(Pic_name,0,0,1);
