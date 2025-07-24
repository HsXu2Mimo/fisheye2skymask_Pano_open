
function cluster_max_p = max_sat_each(max_need_vector, cluster_id)
    % Find unique cluster IDs and their indices
    [C, ~, IC] = unique(cluster_id, 'stable');
    
    % Initialize the output matrix
    cluster_max_p = nan(length(C), 3);
    
    % Preallocate vectors for storing max values and indices
    max_vals = zeros(size(C));
    max_indices = zeros(size(C));
    
    % Loop over each unique cluster ID
    for idk = 1:length(C)
        % Find indices where cluster ID matches the current unique ID
        av_f = IC == idk;
        
        % Get the corresponding values from max_need_vector
        cc_vec = max_need_vector(av_f);
        
        % Find the maximum value and its index in the current cluster
        [M, I] = max(cc_vec);
        
        % Store the results in the preallocated vectors
        max_vals(idk) = M;
        max_indices(idk) = find(av_f, 1, 'first'); % Find the first occurrence of the maximum value
    end
    
    % Populate the output matrix
    cluster_max_p(:, 1) = C;
    cluster_max_p(:, 2) = max_vals;
    cluster_max_p(:, 3) = max_indices;
end













% function cluster_max_p=max_sat_each(max_need_vector, cluster_id)
% max_need_vector=nlos_set; cluster_id=nlos_cluster;
% for idx=
% ismember()
% arrayfun(@(x) mean(max_need_vector(cluster_id==x,1:2),1) ,unique(cluster_id) )
% [C,~,IC] = unique(cluster_id,'rows','stable');
% unimax_need_vector=max_need_vector(IA,:);
% [n,~,~] = histcounts(IC,);%
% REPAET_num=find(n>1);
% Unique_num=find(n==1);%
% cluster_max_p= cell2mat(reshape(arrayfun(@(x) , ,'UniformOutput' ,false),[],1));
% unimax_need_vector(IC(REPAET_num),:)=[];
% cluster_max_p= [cluster_max_p;unimax_need_vector];
% cluster_max_p=cluster_max_p(~isnan(cluster_max_p(:,1)),:);
% cluster_max_p=sortrows(cluster_max_p,1    );
% cluster_max_p=[C,cluster_max_p];
% max(ans)
% REPAET_num=1:max(IC);cluster_max_p=[];
% for idk=1:length(REPAET_num)
%     [M,I]= max(max_need_vector(cluster_id==C( REPAET_num(idk)),:));
%     av_f=find(cluster_id==C( REPAET_num(idk)));
%     cluster_max_p=[cluster_max_p;[C(idk),M,av_f(I)]];
% end
% end




function cluster_max_p=max_sat_each0(max_need_vector, cluster_id)

[C,~,IC] = unique(cluster_id,'rows','stable');

idx_Line_no=1:length(cluster_id);
% par
REPAET_num=1:max(IC);cluster_max_p=[];

% [n,~,~] =histcounts(IC ,'Normalization','count'....
%     ,'BinWidth',1,'BinMethod','integers');
%
% all_ele_tab=nan(length(C),max(n) );
% % -55*ones
% all_ele_tab_idx=nan(length(C),max(n) );
%
% for idk=1:length(C)
% cc_idx=(IC==C);%idx_Line_no
% all_ele_tab(1:sum(cc_idx),idk)=max_need_vector(cc_idx);
%
% all_ele_tab_idx(1:sum(cc_idx),idk)=idx_Line_no(cc_idx);
% end
%
% [M,I]= max(all_ele_tab,[],1,'omitnan');
cluster_max_p=nan(length(REPAET_num),3);


% di

for idk=1:length(REPAET_num)

    %     % av_f=find(cluster_id==C( REPAET_num(idk)));

    %     %     re_cc=C( REPAET_num(idk));
    %     %     av_f=(cluster_id==re_cc);

    av_f= IC==idk;
    %     av_f=logical(IC==idk) ;

    % av_f=ismember_mex( IC, idk );
    %     [av_f] = ismember(IC,idk);
    % av_f=IC==idk;
    cc_line=idx_Line_no(av_f);
    % cc_line=find(av_f);
    cc_vec=max_need_vector(av_f,1);



    [M,I]= max(cc_vec);

    %     cluster_max_p=[cluster_max_p;C(idk),M,cc_line(I)];
    cluster_max_p(idk,1:3)=[C(idk),M,cc_line(I)];

    %     max_need_vector_cc=max_need_vector;
    % %  max_need_vector_cc(IC==idk)
    %     max_need_vector_cc(IC~=idk) =nan;
    %     [M,I]= max(max_need_vector_cc,[],1,'omitnan');
    %     cluster_max_p=[cluster_max_p;[C(idk),M,(I)]];
    %
end


end




















function cluster_max_p=max_sat_each_old(max_need_vector, cluster_id)

[C,~,IC] = unique(cluster_id,'rows','stable');

idx_Line_no=1:length(cluster_id);
% par
REPAET_num=1:max(IC);cluster_max_p=[];
for idk=1:length(REPAET_num)
    
    % av_f=find(cluster_id==C( REPAET_num(idk)));
    
    re_cc=C( REPAET_num(idk));
    av_f=(cluster_id==re_cc);
    
    cc_line=idx_Line_no(av_f);
    
    cc_vec=max_need_vector(av_f,:);
    [M,I]= max(cc_vec);
    
    
    cluster_max_p=[cluster_max_p;[C(idk),M,cc_line(I)]];
end


end