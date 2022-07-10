% ###### SRO_Demo#######

%% Format of the input file
%LAMMPS data file via write_data, version 29 Oct 2020, timestep = 1000067

% 4000 atoms
% 3 atom types
% 
% -0.21379752262759766 35.77379752262689 xlo xhi
% -0.22310728305595262 35.78310728306221 ylo yhi
% -0.1510382480924335 35.71103824808776 zlo zhi
% 
% Masses
% 
% 1 58.69
% 2 58.93
% 3 52
% 
% Atoms # atomic
% 
% 3304 3 1.4643398829714263 4.0928195196098764 3.6047160076132676 -1 -2 2
% 3264 3 1.4547145286449408 0.4285944782019345 3.6111714888843913 -1 -2 2
% 3303 1 3.2949407691122214 2.2796560535324972 3.638080136709751 -1 -2 2

clc;
clear;
addpath(genpath([pwd,'/SRO_module']));
disp(" To show how to use the function to calculate the SRO for a vcmc data")
disp("input: the vcmc output file(the format please go to the help file) ")


lattice_constant=3.556;
supersize=10;
nn_th_SRO=1;
%% read data
data_file_path=['/home/biaoxu4/inflence_of_SRO/structure/950k/'];
% data_file_path=[pwd,'\550k\data.vcmc550'];
% data_file_path=[pwd,'\350k\data.vcmc350'];
[data,infor_set] = read_lmp([data_file_path,'/data.vcmc950']);
index=data(1,:)';
types=data(2,:)';
%random types
% load('random_NiCrCo_types.mat','alloy_type');
% % types=alloy_type';
% 
% [data1,infor_set1]= read_lmp('NiCoCr.lmp');
% data1(2,:)=alloy_type;
% write_lmp('NiCrCo_random',data1',infor_set1)
coords=data(3:5,:)';
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
disp(SRO_set)
  


