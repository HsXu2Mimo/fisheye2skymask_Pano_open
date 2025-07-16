function J_OUT_up=img_Lout_clear(Half_J_OUT,Lout,ex_color)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

if ~isempty(ex_color)
stg_color=ex_color;
clor_resloution=30*1024;
[row_s,col_s,~] = find(Lout~=stg_color);
backgroundInd= sub2ind(size(Lout),row_s,col_s);
[Alab,map_obs] = rgb2ind(Half_J_OUT,clor_resloution);

perm_n=cumprod(size(Alab)   );

perm_n=1:perm_n(2);
perm_n(backgroundInd)=[];

Alab(perm_n)=nan;

J_OUT_up=ind2rgb(Alab,map_obs);
else
  J_OUT_up=  Half_J_OUT;
end
end