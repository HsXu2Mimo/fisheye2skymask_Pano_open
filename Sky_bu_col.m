function   [master_color,cand_color,top_color]=Sky_bu_col(Lout)

% Copyright (C) 2020-2024 Xuhaosheng
% All rights reserved.
% 17074845g@connect.polyu.hk

[L_m,L_n]=size(Lout);
numColors=length(unique(Lout));
cand_color=[];

cand_color(:,1)=[unique(Lout)];
for  co_id=1:numColors

    [i,j]=find(Lout==cand_color(co_id,1));

    cand_color(co_id,[2:3,4])=[mean([i,j],1),(sum(Lout==cand_color(co_id,1),'all')./L_m./L_n)];

end
[~,I]=max(cand_color(:,2));
master_color=cand_color(I,1);

counter_bu_1=(sum(Lout==master_color,'all')./L_m./L_n);
[~,I]=min(cand_color(:,2));
top_color=cand_color(I,1);

counter_bu_2=(sum(Lout==top_color,'all')./L_m./L_n);


counter_bu=(sum(Lout~=master_color,'all')./L_m./L_n);
disp(['?!!Building Inside Tot->-' , '->Color-Mode->' num2str(counter_bu)])

disp(['?!!Building Down ->-' ,num2str(master_color) '/P--' num2str(counter_bu_1)])
disp(['?!!Building tOP ->-' , num2str(top_color) '/P--' num2str(counter_bu_2)])

cand_color=sortrows(cand_color,2);
end