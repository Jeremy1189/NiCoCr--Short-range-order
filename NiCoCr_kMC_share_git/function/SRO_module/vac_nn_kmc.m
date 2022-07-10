function [vac_ID_nn_set,vac_ID_count_set]= vac_nn_kmc(vac_ID,per_coords,lattice_constant,supersize,largest_nn)
% Description:
%       This function aims to obtain the atomIDs in nn shells
% input:
%      vac_ID:size=[1,1]
%      mig_ID: size=[1,1]
%      per_coords: size=[n,3]
%       largest_nn= integer,1,2,...,n
% output:
%       vac_mig_sortID_nn_set: 
%       vac_mig_nn_count_set
%       vac_ID_nn_set:
%       vac_ID_count_set:
% time:
%          2020/11/15 (first version)
%%
vac_coord = per_coords(vac_ID,:);
L = length(per_coords);
relative_vac_coords = per_coords - repmat(vac_coord,[L,1]);
fix_point_trans=10^8;
relative_vac_coords= round(relative_vac_coords.*fix_point_trans)./fix_point_trans;
% periodic boundary check
min_per_coords= min(per_coords);
% max_per_coords= max(per_coords);
% lattice_constant=3.488;
% supersize=10;
max_per_coords=lattice_constant*supersize*ones(1,3);
box_size =  max_per_coords-min_per_coords;
upper_boundary = round(1/2.*(box_size+lattice_constant/2).*fix_point_trans)./fix_point_trans;

% err_tol =10^-2;
for i= 1:length(box_size)
    vac_index_modify= find(relative_vac_coords(:,i)>upper_boundary(i));
   
%     vac_index_modify= find(relative_vac_coords(:,i)-1/2*box_size(i)>err_tol);
%     mig_index_modify= find( relative_mig_coords(:,i)-1/2*box_size(i)>err_tol);

    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)-box_size(i);
    % less than and equal to the box_size
%     vac_index_modify= find(relative_vac_coords(:,i)+1/2*box_size(i)<-err_tol);
%     mig_index_modify= find( relative_mig_coords(:,i)+1/2*box_size(i)<-err_tol);
    vac_index_modify= find(relative_vac_coords(:,i)<=-upper_boundary(i));
    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)+box_size(i);
end
% calculated the distance
vac_nn_distance = sqrt( sum(relative_vac_coords.^2,2) );
% fixed
fixed_point=8;
vac_nn_distance = round(vac_nn_distance*10^fixed_point)/10^fixed_point;
%unique distance
%vac
unique_vac_dis_set= unique(vac_nn_distance);
% largest_nn=6;
vac_ID_nn_set =zeros(largest_nn,1);% output1
vac_ID_count_set= zeros(largest_nn,1);% remove vacancy itself distance 
% L_unique = length(unique_mig_dis_set); % it is euqal with the above L_unique

count_vac=0;
for num_unique = 1:largest_nn
    %vac
    index_vac_set=find(vac_nn_distance ==unique_vac_dis_set(num_unique+1));% the first is itself
    vac_ID_count_set(num_unique) = length(index_vac_set);
    vac_ID_nn_set(count_vac+1 : count_vac+ length(index_vac_set) ) = index_vac_set;%output2
    count_vac = count_vac + length(index_vac_set);
end
% output

end

