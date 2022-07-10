function Central_ID_nn_set = vac_coords_trans(central_coord,per_coords,per_index,lattice_constant,supersize,nn_set)
% Description:
%       This function aims to obtain the atomIDs in nn shells
% input:
%      vac_coord:size=[1,3]
%       per_coords: size=[n,3]

% output:
%       vac_ID_nn_set:
%       vac_ID_count_set:
% time:
%          2020/11/11 (first version)
%%
relative_vac_coords = per_coords - repmat(central_coord,[length(per_coords),1]);
fix_point_trans=10^8;
relative_vac_coords= round(relative_vac_coords.*fix_point_trans)./fix_point_trans;
% periodic boundary check
min_per_coords= zeros(1,3);
% lattice_constant=3.522;
% supersize=10;
max_per_coords=lattice_constant*supersize*ones(1,3);
box_size =  max_per_coords-min_per_coords;
upper_boundary = round(1/2.*box_size.*fix_point_trans)./fix_point_trans+lattice_constant/4;
% err_tol =10^-2;
for i= 1:length(box_size)
    vac_index_modify= find(relative_vac_coords(:,i)>upper_boundary(i));

    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)-box_size(i);

    vac_index_modify= find(relative_vac_coords(:,i)<=-upper_boundary(i));
  
    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)+box_size(i);
   
end
% calculated the distance
vac_nn_distance = sqrt( sum(relative_vac_coords.^2,2) );

% fixed
fixed_point=8;
vac_nn_distance = round(vac_nn_distance*10^fixed_point)/10^fixed_point;
% 
[nn_values,nn_index]=sort(vac_nn_distance,'ascend');
nn_index(1)=[];% remove itself
sort_per_index=per_index(nn_index);

cum_nn_set=cumsum(nn_set);
Central_ID_nn_set=[];
for j=1:length(cum_nn_set)
    if j==1
        nn_count=[1:cum_nn_set(j)]';
    else
        nn_count=[1+cum_nn_set(j-1):cum_nn_set(j)]';
    end
    Central_ID_nn_set=[Central_ID_nn_set;sort_per_index(nn_count)];
    
end

%

end%end function

