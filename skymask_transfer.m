function skymask=skymask_transfer(skymask,heading_rotation,mirror_1)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

if nargin<=2
    mirror_1=0;
end
if mirror_1==1
% % mirror
skymask(:,1)=-skymask(:,1)+360;
skymask=sortrows(skymask,1);
end
skymask_1_save=[];
skymask_1_save(:,1)=skymask(:,1);
% % rotation

heading_rotation=mod(heading_rotation+360*2,360);
if isempty(heading_rotation)
heading_rotation=0;
end
if isnan(heading_rotation)

   skymask=nan(361,2);
   return
end
skymask(:,1)=skymask(:,1)+heading_rotation;

skymask((skymask(:,1)>360),1)=skymask((skymask(:,1)>360),1)-360;
skymask=sortrows(skymask,1);
% skymask(:,1)=mod(skymask(:,1),360);
% skymask=sortrows(skymask,1);
sky_temp=....
    mean_sat_each(skymask(:,2),skymask(:,1));
if ~isnan(sky_temp(:,1)) 

end
f0=fit(sky_temp(~isnan(sky_temp(:,2)) ,1),sky_temp(~isnan(sky_temp(:,2)) ,2),'linearinterp');


skymask_1_save(~isnan(sky_temp(:,2)),2)=feval(f0,skymask_1_save(~isnan(sky_temp(:,2)),1));

skymask_1_save(isnan(sky_temp(:,2)),2)=nan;

skymask=skymask_1_save;

if  skymask(skymask(:,1)==0,2)~=skymask(skymask(:,1)==360,2)
% clc

tp_el=max([skymask(skymask(:,1)==0,2),skymask(skymask(:,1)==360,2)]);

skymask(logical((skymask(:,1)==0)+....
    (skymask(:,1)==360))      ,2) =tp_el;

end
end
function cluster_max_p=mean_sat_each(max_need_vector, cluster_id)

[C,~,IC] = unique(cluster_id,'rows','stable');

idx_Line_no=1:length(cluster_id);
% par
REPAET_num=1:max(IC);cluster_max_p=[];

% cluster_max_p=nan(length(REPAET_num),3);
for idk=1:length(REPAET_num)
    av_f= IC==idk;
    cc_line=idx_Line_no(av_f);
    % cc_line=find(av_f);
    cc_vec=max_need_vector(av_f,1);
    [M]= mean(cc_vec,1);

    cluster_max_p(idk,:)=[C(idk),M,cc_line(1)];
    %
end

end
