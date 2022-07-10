% function  SRO_table=SRO_cal_NiCrCo(crystal_structure_n_5,statist_SRO_nn_1_1)
function  SRO_set=SRO_cal_NiCrCo(coords,types,lattice_constant,supersize,nn_th_SRO)
%input: 
% crystal_structure:[atom_index, atom_types, atom_coordinate_x,  atom_coordinate_y, atom_coordinate_z]
% date:2022/1/23
% output:
%  SRO_table= [key_types: SRO values]
% Algorithm reference: Warren-Cowly SRO alpah=1-p(x,T)/x , x is the composition, T is the temperature. p(x,T) is the
% condition probability that given an A atom at the origin, there is a B


num_total_atoms=length(types);
coordnate_nn=[12,6,24,12,24,8];
coordinate_num_atoms=coordnate_nn(nn_th_SRO)*num_total_atoms;


unique_types=unique(types);
L_types =length(unique_types);
c_set=zeros(L_types,1);
for j=1:L_types
    cur_index=find(types==unique_types(j));
    c_set(j)= length(cur_index)/num_total_atoms;   
end

%% category
NN_type_cell=cell(L_types,1);
SRO_set=zeros(L_types,L_types);
for num_types=1:L_types
    cur_index=find(types==unique_types(num_types));
    L_cur=length(cur_index);
    cur_type_coords=coords(cur_index,:);
    % the center NN1 of the current 
    cur_nn_type_set=zeros( L_cur,200);% just give a large inital matrix
    for num_type_count=1:L_cur
        central_coord=cur_type_coords(num_type_count,:);
        [dis_value,index_nn] = nn_dis_count(central_coord,coords,lattice_constant,supersize,nn_th_SRO);
      
        for s=1:length(index_nn)
            cur_nn_type_set(num_type_count,s)=types(index_nn(s));
        end
    end
    NN_type_cell{num_types} =  cur_nn_type_set;
    c_cur= L_cur/num_total_atoms;
    for k=1:L_types
        cur_pairs_index= find(cur_nn_type_set==unique_types(k));
        num_cur_pairs=length(cur_pairs_index);
        prob_cur_pair=num_cur_pairs/(coordinate_num_atoms*c_cur);
        
        SRO_set(num_types,k)=1-prob_cur_pair/c_set(k);
    end
    
    
end
% disp(SRO_set)