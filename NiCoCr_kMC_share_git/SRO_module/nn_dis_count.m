function [dis_value,index_nn] = nn_dis_count(central_coord,per_coords,lattice_constant,supersize,nn_th_SRO)
% Derroescription:
%       This function aims to obtain the atomIDs in nn shells
% input:


% output:

% time:
%          2022/01/23 (first version)

%% obtain the cutoff value of the n_th NN
if nn_th_SRO<0 || mod(nn_th_SRO,1)~=0
    error('nn_th_SRO must be integer!!!!')
    
end
unit=sqrt(2)/2*lattice_constant;%fcc 1nn distance
 nn_distanse= zeros(nn_th_SRO,1);
for i=1:nn_th_SRO+1
    nn_distanse(i)= sqrt(i)*unit;       
end
% cut_off value is the less than (nn_distance+half_distantce_between_nn_and
% nn+1), and it large than (nn_distance-half_distantce_between_nn_and nn-1)
half_dis_between_nn_layers=(nn_distanse(2:end)-nn_distanse(1:end-1))./2;
cut_off_dis_upper=nn_distanse(end-1)+half_dis_between_nn_layers(end);% nn and nn+1
if nn_th_SRO==1
    cut_off_dis_lower=nn_distanse(end-1)-nn_distanse(end-1)/2;
else
    cut_off_dis_lower=nn_distanse(end-1)-half_dis_between_nn_layers(end-1);%nn and nn-1
end

%% relative coodinates transfer    
relative_vac_coords = per_coords - repmat(central_coord,[length(per_coords),1]);
fix_point_trans=10^8;
relative_vac_coords= round(relative_vac_coords.*fix_point_trans)./fix_point_trans;
% periodic boundary check,
min_per_coords= zeros(1,3);
% lattice_constant=3.522;
% supersize=10;
max_per_coords=lattice_constant*supersize*ones(1,3);
box_size =  max_per_coords-min_per_coords;
upper_boundary_each_dim = round(1/2.*box_size.*fix_point_trans)./fix_point_trans+lattice_constant/4;
% err_tol =10^-2;
for i= 1:length(box_size)
    vac_index_modify= find(relative_vac_coords(:,i)>upper_boundary_each_dim(i));

    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)-box_size(i);

    vac_index_modify= find(relative_vac_coords(:,i)<=-upper_boundary_each_dim(i));
  
    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)+box_size(i);
   
end
% calculated the distance
vac_nn_distance = sqrt( sum(relative_vac_coords.^2,2) );

% fixed
fixed_point=8;
vac_nn_distance = round(vac_nn_distance*10^fixed_point)/10^fixed_point;
% 
% [nn_values,nn_index]=sort(vac_nn_distance,'ascend');

   index_nn = find(vac_nn_distance>cut_off_dis_lower & vac_nn_distance<cut_off_dis_upper);

   dis_value=vac_nn_distance(index_nn);

end



