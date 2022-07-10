function  SRO_table=SRO_cal_NiCrCo(crystal_structure,vac_ID,supersize,lattice_constant,nn)
%input: 
% crystal_structure:[atom_index, atom_types, atom_coordinate_x,  atom_coordinate_y, atom_coordinate_z]
% date:2021/10/29
% output:
%  SRO_table= [key_types: SRO values]
% Algorithm reference: Warren-Cowly SRO alpah=1-p(x,T)/x , x is the composition, T is the temperature. p(x,T) is the
% condition probability that given an A atom at the origin, there is a B
% atom at shell
index = crystal_structure(:,1);
types = crystal_structure(:,2);
coords=crystal_structure(:,3:5);



%% category
Ni_index=find(types==1);
Co_index=find(types==2);
Cr_index=find(types==3);
Ni_coords= coords(Ni_index,:); %#ok<*FNDSB>
Co_coords= coords(Co_index,:);
Cr_coords= coords(Cr_index,:);
%% neighbours for each Ni atom
% load neigh_count_set.mat neigh_count_set
% nn=5;
% vac_ID=2221;
% supersize=10;
% lattice_constant=3.556;
[~,vac_ID_count_set]= vac_nn_kmc(vac_ID,coords,lattice_constant,supersize,nn+2);
neigh_count_set=vac_ID_count_set;
nn_set=neigh_count_set(1:nn);
cum_nn_set=cumsum(nn_set);
%% Ni atom as the central

if nn==1
    s_p=1;
elseif mod(nn,1)==0
    s_p=cum_nn_set(nn-1)+1;
else
    error('invalid input nn value, it must be a integer');
end
e_p=cum_nn_set(nn);
Ni_nn_type_set=zeros(e_p-s_p+1,length(Ni_coords));
for num_Ni=1:length(Ni_coords)
    central_coord=Ni_coords(num_Ni,:);
    Central_ID_nn_set= vac_coords_trans_sro(central_coord,coords,index,lattice_constant,supersize,nn_set);
    for s=1:(e_p-s_p+1)
        index_s=find(index==Central_ID_nn_set(s),1);
        Ni_nn_type_set(s,num_Ni)=types(index_s); 
    end
           
end
%% Co atom as the cetral
Co_nn_type_set=zeros(e_p-s_p+1,length(Ni_coords));
for num_Co=1:length(Co_coords)
    central_coord=Co_coords(num_Co,:);
    Central_ID_nn_set= vac_coords_trans_sro(central_coord,coords,index,lattice_constant,supersize,nn_set);
    for s=1:(e_p-s_p+1)
        index_s=find(index==Central_ID_nn_set(s),1);
        Co_nn_type_set(s,num_Co)=types(index_s); 
    end
           
end
%%
%% Cr atom as the cetral
Cr_nn_type_set=zeros(e_p-s_p+1,length(Ni_coords));
for num_Cr=1:length(Cr_coords)
    central_coord=Cr_coords(num_Cr,:);
    Central_ID_nn_set= vac_coords_trans_sro(central_coord,coords,index,lattice_constant,supersize,nn_set);
    for s=1:(e_p-s_p+1)
        index_s=find(index==Central_ID_nn_set(s),1);
        Cr_nn_type_set(s,num_Cr)=types(index_s); 
    end
           
end
%% calculate SRO with central atom Ni
% c_Ni=1/3;
% c_Co=1/3;
% c_Cr=1/3;
c_Ni=length(Ni_index)/length(index);
c_Co=length(Co_index)/length(index);
c_Cr=length(Cr_index)/length(index);
num_total_atoms=length(types);


Ni_NN1_types_set= Ni_nn_type_set;
total_nn=size(Ni_NN1_types_set,1)*num_total_atoms;
%Ni-Ni
Ni_Ni_index=find(Ni_NN1_types_set==1);
num_Ni_Ni=length(Ni_Ni_index);% each Ni-Ni pairs counted two times
P_Ni_Ni=num_Ni_Ni/(total_nn);% probability of Ni-Ni pairs
SRO_Ni_Ni_NN1=1- P_Ni_Ni/(c_Ni*c_Ni);
% Ni-Co
Ni_Co_index=find(Ni_NN1_types_set==2);
num_Ni_Co=length(Ni_Co_index);% each Ni-Ni pairs counted two times
% P_Ni_Co=num_Ni_Co/(total_nn/2);% probability of Ni-Co pairs, please go to the reference ppt for the formular
P_Ni_Co=num_Ni_Co/(total_nn);
SRO_Ni_Co_NN1=1- P_Ni_Co/(c_Ni*c_Co);
% Ni-Cr
Ni_Cr_index=find(Ni_NN1_types_set==3);
num_Ni_Cr=length(Ni_Cr_index);% each Ni-Ni pairs counted two times
% P_Ni_Cr=num_Ni_Cr/(total_nn/2);% probability of Ni-Co pairs
P_Ni_Cr=num_Ni_Cr/(total_nn);
SRO_Ni_Cr_NN1=1- P_Ni_Cr/(c_Ni*c_Cr);
%% calculate SRO with central atom Co
Co_NN1_types_set= Co_nn_type_set;
%Co-Co
Co_Co_index=find(Co_NN1_types_set==2);
num_Co_Co=length(Co_Co_index);% 
P_Co_Co=num_Co_Co/(total_nn);% probability of Co-Co pairs
SRO_Co_Co_NN1=1- P_Co_Co/(c_Co*c_Co);
% Co-Ni
Co_Ni_index=find(Co_NN1_types_set==1);
num_Co_Ni=length(Co_Ni_index);% each Ni-Ni pairs counted two times
P_Co_Ni=num_Co_Ni/(total_nn);% probability of Ni-Co pairs
SRO_Co_Ni_NN1=1- P_Co_Ni/(c_Co*c_Ni);
% Co-Cr
Co_Cr_index=find(Co_NN1_types_set==3);
num_Co_Cr=length(Co_Cr_index);% each Ni-Ni pairs counted two times
P_Co_Cr=num_Co_Cr/(total_nn);% probability of Ni-Co pairs
SRO_Co_Cr_NN1=1- P_Co_Cr/(c_Co*c_Cr);
%% calculate SRO with central atom Cr
Cr_NN1_types_set= Cr_nn_type_set;
%Cr-Cr
Cr_Cr_index=find(Cr_NN1_types_set==3);
num_Cr_Cr=length(Cr_Cr_index);% each Ni-Ni pairs counted two times
P_Cr_Cr=num_Cr_Cr/(total_nn);% probability of Ni-Ni pairs
SRO_Cr_Cr_NN1=1- P_Cr_Cr/(c_Cr*c_Cr);
% Cr-Ni
Cr_Ni_index=find(Cr_NN1_types_set==1);
num_Cr_Ni=length(Cr_Ni_index);% each Ni-Ni pairs counted two times
P_Cr_Ni=num_Cr_Ni/(total_nn);% probability of Ni-Co pairs
SRO_Cr_Ni_NN1=1- P_Cr_Ni/(c_Cr*c_Ni);
% Cr-Co
Cr_Co_index=find(Cr_NN1_types_set==2);
num_Cr_Co=length(Cr_Co_index);% each Ni-Ni pairs counted two times
P_Cr_Co=num_Cr_Co/(total_nn);% probability of Ni-Co pairs
SRO_Cr_Co_NN1=1- P_Cr_Co/(c_Cr*c_Co);

T= table(SRO_Ni_Ni_NN1,SRO_Ni_Co_NN1,SRO_Ni_Cr_NN1,...
    SRO_Co_Co_NN1,SRO_Co_Ni_NN1,SRO_Co_Cr_NN1,...
    SRO_Cr_Cr_NN1,SRO_Cr_Ni_NN1,SRO_Cr_Co_NN1);

%% output
SRO_table=T;