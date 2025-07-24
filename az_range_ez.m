

function [near_sat, pick_idex] = az_range_ez(azi_range, near_sat_time_prn_az_ele, az_id)

    if nargin < 3
        az_id = 3;
    end

    azi_range = sort(azi_range);
    azi_range = reshape(azi_range, 1, []);

    if length(azi_range) == 1
        azi_range = [azi_range - 1, azi_range + 1];
    end

    if nargin > 2 && numel(near_sat_time_prn_az_ele) == 1
        near_sat_time_prn_az_ele = reshape(near_sat_time_prn_az_ele, [], 1);
        near_sat_time_prn_az_ele = repmat(near_sat_time_prn_az_ele, 1, length(near_sat_time_prn_az_ele));
    end

    near_sat = near_sat_time_prn_az_ele;

    if any(azi_range < 0)
        pick_idex_1 = find(near_sat(:, az_id) >= azi_range(1) + 360 & near_sat(:, az_id) < 360);
        pick_idex_2 = find(near_sat(:, az_id) >= 0 & near_sat(:, az_id) <= azi_range(2));
        pick_idex = [pick_idex_1; pick_idex_2];
    else
        pick_idex = find(near_sat(:, az_id) >= azi_range(1) & near_sat(:, az_id) <= azi_range(2));
    end

    near_sat = near_sat(pick_idex, :);
end











function [near_sat,pick_idex]=az_range_ez0(azi_range,near_sat_time_prn_az_ele,az_id)

if nargin<=2
    az_id=3;
    %     mod,360
end
azi_range=sort( (azi_range));
azi_range=reshape( sort(azi_range),1,[])   ;
if length(azi_range)==1
    azi_range=[azi_range-1,azi_range+1];
    %     disp(['azi_range' num2str(azi_range)])
end


if nargin>2&& sum(size( near_sat_time_prn_az_ele )==[1])>0
    %     clc
    near_sat_time_prn_az_ele=reshape(near_sat_time_prn_az_ele,[],1);
    near_sat_time_prn_az_ele=repmat(near_sat_time_prn_az_ele,1,length(near_sat_time_prn_az_ele));
    %      disp(['azi_range ---->size ' num2str(size( near_sat_time_prn_az_ele ))])
end


near_sat=near_sat_time_prn_az_ele;
pick_idex=find(logical([near_sat(:,az_id)>=azi_range(1)].*......
    [near_sat(:,az_id)<=azi_range(2)]));

if any(azi_range<0)
    pick_idex_1=find(logical([near_sat(:,az_id)>=azi_range(1)+360].*......
        [near_sat(:,az_id)<360]));
    pick_idex_2=find(logical([near_sat(:,az_id)>=0].*......
        [near_sat(:,az_id)<=azi_range(2)]));
    %     disp(['azi_range' num2str(azi_range),'-- Less 0'])
    %     disp(['New_azi_range' num2str(azi_range(1)+360),'-','/360/' num2str(num2str(azi_range(2)))])
    pick_idex=[pick_idex_1;pick_idex_2];
else
    %    disp(['azi_range' num2str(azi_range),'-- Less 0'])
    %    disp(['New_azi_range' num2str(azi_range(1)+360),'-',''])
    pick_idex=find(logical([near_sat(:,az_id)>=azi_range(1)].*......
        [near_sat(:,az_id)<=azi_range(2)]));
end



% end

% sum(logical([near_sat(:,az_id)>=azi_range(1)].*......
%     [near_sat(:,az_id)<=azi_range(2)]))
% near_sat=near_sat( logical([near_sat(:,az_id)>=azi_range(1)].*......
%     [near_sat(:,az_id)<=azi_range(2)]),:);
% near_sat=sortrows(near_sat,az_id);

near_sat=near_sat(pick_idex,:);




end