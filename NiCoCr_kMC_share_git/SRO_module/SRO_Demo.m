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

%% category
Ni_index=find(types==1);
Co_index=find(types==2);
Cr_index=find(types==3);
% Ni_types=types(Ni_index);
% Co_types=types(Co_index);
% Cr_types=types(Cr_index);
Ni_coords= coords(Ni_index,:);
Co_coords= coords(Co_index,:);
Cr_coords= coords(Cr_index,:);
% Ni atoms
Ni_nn_type_set=zeros(length(Ni_coords),200);% just give a large inital matrix
count=0;
for num_Ni=1:length(Ni_coords)
    central_coord=Ni_coords(num_Ni,:);
    [dis_value,index_nn] = nn_dis_count(central_coord,coords,lattice_constant,supersize,nn_th_SRO);


    for s=1:length(index_nn)

        Ni_nn_type_set(num_Ni,s)=types(index_nn(s));
    end
           
end

Co_nn_type_set=zeros(length(Co_coords),200);% just give a large inital matrix

for num_Co=1:length(Co_coords)
    central_coord=Co_coords(num_Co,:);
    [dis_value,index_nn] = nn_dis_count(central_coord,coords,lattice_constant,supersize,nn_th_SRO);

    for s=1:length(index_nn)

        Co_nn_type_set(num_Co,s)=types(index_nn(s));
    end
           
end

Cr_nn_type_set=zeros(length(Cr_coords),200);% just give a large inital matrix
for num_Cr=1:length(Cr_coords)
    central_coord=Cr_coords(num_Cr,:);
    [dis_value,index_nn] = nn_dis_count(central_coord,coords,lattice_constant,supersize,nn_th_SRO);

    for s=1:length(index_nn)

        Cr_nn_type_set(num_Cr,s)=types(index_nn(s));
    end
           
end
c_Ni=length(Ni_index)/length(index);
c_Co=length(Co_index)/length(index);
c_Cr=length(Cr_index)/length(index);
num_total_atoms=length(types);
Ni_NN1_types_set= Ni_nn_type_set;
total_nn1_N=12*num_total_atoms;
%Ni-Ni
Ni_Ni_index=find(Ni_NN1_types_set==1);
num_Ni_Ni=length(Ni_Ni_index);% each Ni-Ni pairs counted two times
P_Ni_Ni=num_Ni_Ni/(total_nn1_N);% probability of Ni-Ni pairs
SRO_Ni_Ni_NN1=1- P_Ni_Ni/(c_Ni*c_Ni);
% Ni-Co
Ni_Co_index=find(Ni_NN1_types_set==2);
num_Ni_Co=length(Ni_Co_index);% each Ni-Ni pairs counted two times
% P_Ni_Co=num_Ni_Co/(total_nn1_N/2);% probability of Ni-Co pairs, please go to the reference ppt for the formular
P_Ni_Co=num_Ni_Co/(total_nn1_N);
SRO_Ni_Co_NN1=1- P_Ni_Co/(c_Ni*c_Co);
% Ni-Cr
Ni_Cr_index=find(Ni_NN1_types_set==3);
num_Ni_Cr=length(Ni_Cr_index);% each Ni-Ni pairs counted two times
% P_Ni_Cr=num_Ni_Cr/(total_nn1_N/2);% probability of Ni-Co pairs
P_Ni_Cr=num_Ni_Cr/(total_nn1_N);
SRO_Ni_Cr_NN1=1- P_Ni_Cr/(c_Ni*c_Cr);
%% calculate SRO with central atom Co
Co_NN1_types_set= Co_nn_type_set;
%Co-Co
Co_Co_index=find(Co_NN1_types_set==2);
num_Co_Co=length(Co_Co_index);% each Ni-Ni pairs counted two times
P_Co_Co=num_Co_Co/(total_nn1_N);% probability of Ni-Ni pairs
SRO_Co_Co_NN1=1- P_Co_Co/(c_Co*c_Co);
% Co-Ni
Co_Ni_index=find(Co_NN1_types_set==1);
num_Co_Ni=length(Co_Ni_index);% each Ni-Ni pairs counted two times
P_Co_Ni=num_Co_Ni/(total_nn1_N);% probability of Ni-Co pairs
SRO_Co_Ni_NN1=1- P_Co_Ni/(c_Co*c_Ni);
% Co-Cr
Co_Cr_index=find(Co_NN1_types_set==3);
num_Co_Cr=length(Co_Cr_index);% each Ni-Ni pairs counted two times
P_Co_Cr=num_Co_Cr/(total_nn1_N);% probability of Ni-Co pairs
SRO_Co_Cr_NN1=1- P_Co_Cr/(c_Co*c_Cr);
%% calculate SRO with central atom Cr
Cr_NN1_types_set= Cr_nn_type_set;
%Cr-Cr
Cr_Cr_index=find(Cr_NN1_types_set==3);
num_Cr_Cr=length(Cr_Cr_index);% each Ni-Ni pairs counted two times
P_Cr_Cr=num_Cr_Cr/(total_nn1_N);% probability of Ni-Ni pairs
SRO_Cr_Cr_NN1=1- P_Cr_Cr/(c_Cr*c_Cr);
% Cr-Ni
Cr_Ni_index=find(Cr_NN1_types_set==1);
num_Cr_Ni=length(Cr_Ni_index);% each Ni-Ni pairs counted two times
P_Cr_Ni=num_Cr_Ni/(total_nn1_N);% probability of Ni-Co pairs
SRO_Cr_Ni_NN1=1- P_Cr_Ni/(c_Cr*c_Ni);
% Cr-Co
Cr_Co_index=find(Cr_NN1_types_set==2);
num_Cr_Co=length(Cr_Co_index);% each Ni-Ni pairs counted two times
P_Cr_Co=num_Cr_Co/(total_nn1_N);% probability of Ni-Co pairs
SRO_Cr_Co_NN1=1- P_Cr_Co/(c_Cr*c_Co);

T= table(SRO_Ni_Ni_NN1,SRO_Ni_Co_NN1,SRO_Ni_Cr_NN1,...
    SRO_Co_Co_NN1,SRO_Co_Ni_NN1,SRO_Co_Cr_NN1,...
    SRO_Cr_Cr_NN1,SRO_Cr_Ni_NN1,SRO_Cr_Co_NN1)