function  SRO_set= SRO_cal(lattice_constant,supersize,nn_th_SRO,data)
%%
% input:
%  lattice_constant: the value of lattice constant of the supercell eg:5.2
%supersize: eg: 10 , whcih means a 10*10*10 supercell contains 4000 atoms in fcc
%nn_th_SRO: 1 means the 1st nearest neigbour(NN) SRO, 2 means  the 2nd NN SRO
%data : the supercell atom information: [index(1*n);types(1*n);Coordinates(3*n)]
% SRO_set:
  % retures a matrix
  % for example: NiCoCr_SRO:(3(Ni,Co,Cr)*3(Ni,Co,Cr)), which means the SRO_set(1,1)--NI-Ni SRO_set(2,1)--Ni-Co SRO_set(3,1)--Ni-Cr
  %                                                                     SRO_set(2,1)--Co-Ni SRO_set(2,2)--Co-Co  ,SRO_set(2,3)--Co-Cr  
  %                                                                     SRO_set(3,1)--Cr-Ni SRO_set(3,2)--Cr-Co  ,SRO_set3,3)--Cr-Cr


%%
types=data(2,:)';
coords=data(3:5,:)';
num_total_atoms=length(types);
coordnate_nn=[12,6,24,12,24,8];
% cum_coordnate_nn=cumsum(coordnate_nn);
coordinate_num_atoms=coordnate_nn(nn_th_SRO)*num_total_atoms;

%% total  
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
%         for s=1:cum_coordnate_nn(nn_th_SRO)
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