clc;
clear;
close all;

% data_path='/home/biaoxu4/inflence_of_SRO/diffusion/1100k/random_Ni_CoCr_ML_1300k.mat';
data_path='/home/biaoxu4/inflence_of_SRO/diffusion/1100k/restart_Ni_CoCr_ML_1300k.mat';
load(data_path,'R_0', 'R_t', 'ratio', 'lattice_constant_set', 'per_set_ratios', 'all_atom_set_cell', 'tracks_set','per_coords_set','T','iter_times')

%% recover the tracks
Set=[];
supersize=10;
vacancy_atomID=2221;
N=3999;
Lw=2;
D=zeros(length(ratio),1);
D_Co=zeros(length(ratio),1);
D_Ni=zeros(length(ratio),1);
D_Cr=zeros(length(ratio),1);
f_tracer_all=zeros(length(ratio),1);
f_tracer_Co=zeros(length(ratio),1);
f_tracer_Cr=zeros(length(ratio),1);
f_tracer_Ni=zeros(length(ratio),1);
step_sample=1e4;
SRO_cell=cell(length(ratio),round(iter_times/step_sample));
% R_set= cell(length(ratio),round(iter_times/step_sample));
 strcture_cell=cell(length(ratio),round(iter_times/step_sample));
for num_ratio=1:length(ratio)
% for num_ratio=3
    R0=R_0{num_ratio};
    Rt=R_t{num_ratio};
    per=per_set_ratios{num_ratio};
    tracks=tracks_set{num_ratio};
    lattice_constant=lattice_constant_set(num_ratio);
    all_atoms_mig_ID=all_atom_set_cell{num_ratio}(:,1);    
    R_temp=R0;
    cross_boundary_count=zeros(3,length(R0));
    vector_dif_set=zeros(length(tracks)-1,3);
    MSD=zeros(length(tracks),1);
    MSD_Co=zeros(length(tracks),1);
    MSD_Cr=zeros(length(tracks),1);
    MSD_Ni=zeros(length(tracks),1);
    nn=length(tracks)-1; 
    index_set=per(:,1);
    types=per(:,2);
    inital_data=R0;
    inital_data(:,2221)=[];  
    index_set(vacancy_atomID)=[];
    Ni_index=find(types==1);
    Co_index=find(types==2); 
    Cr_index=find(types==3);
    Len=lattice_constant*supersize;
    mig_type_set=zeros(nn,1);
    % for i=1:length(tracks)-1
    
    
     for i=1:nn
%     for i=3
        mig_ID=all_atoms_mig_ID(i);
        mig_type_set(i)=types(mig_ID);
        current_jump_vector= tracks(i,2:4)-tracks(i+1,2:4);
        vector_dif_set(i,:)=  current_jump_vector;
        index_out_L= find(current_jump_vector>=0.5*Len);
        current_jump_vector(index_out_L)=current_jump_vector(index_out_L)-Len;
        index_out_R= find(current_jump_vector<-0.5*Len);
        current_jump_vector(index_out_R)=current_jump_vector(index_out_R)+Len;
        %     disp(['before: mig_ID=',num2str(mig_ID)])
        %     disp([num2str(R_temp(:,mig_ID))])
        next_vac_coord= R_temp(:,mig_ID);
        R_temp(:,mig_ID) =R_temp(:,mig_ID)+current_jump_vector';
        
        %    disp(['after: mig_ID=',num2str(mig_ID)])
        %      disp([num2str(R_temp(:,mig_ID))])
        Rtemp_coords=R_temp(:, mig_ID);
        index_out_U= find(Rtemp_coords>Len);
        cross_boundary_count(index_out_U,mig_ID)= cross_boundary_count(index_out_U,mig_ID)+1;
        Rtemp_coords(index_out_U)= Rtemp_coords(index_out_U)-Len;
        
        index_out_L= find(Rtemp_coords<0);
        cross_boundary_count(index_out_L,mig_ID)= cross_boundary_count(index_out_L,mig_ID)-1;
        Rtemp_coords(index_out_L)= Rtemp_coords(index_out_L)+Len;
        R_temp(:, mig_ID)=Rtemp_coords;
        
        R_r_vacancy = R_temp;
        R_r_vacancy(:,vacancy_atomID)=next_vac_coord;
        if mod(i,step_sample)==0
            count=round(i/step_sample);
            crystal_structure=[per(:,[1,2])';R_r_vacancy];
            strcture_cell{num_ratio,count}=crystal_structure;
            nn_statistic=2;
            SRO_table=SRO_cal_NiCrCo(crystal_structure',vacancy_atomID,supersize,lattice_constant,nn_statistic);
            SRO_cell{num_ratio,count}=SRO_table;
        end
        
        Dt =R_temp-R0;
        Dt=Dt+Len.* cross_boundary_count;       
        sum_Dt=sum(Dt.^2);%
        sum_Dt(vacancy_atomID)=0;
        MSD(i+1) =sum(sum_Dt)./N;
        %partitial MSD
        if ~isempty(Co_index)
            sum_Fe= sum_Dt(Co_index);
            r_index=find(Co_index==vacancy_atomID);
            sum_Fe(r_index)=[];
            MSD_Co(i+1) =sum(sum_Fe)./length(Co_index);
        end
        if ~isempty(Cr_index)
            sum_Fe= sum_Dt(Cr_index);
            r_index=find(Cr_index==vacancy_atomID);
            sum_Fe(r_index)=[];
            MSD_Cr(i+1) =sum(sum_Fe)./length(Cr_index);
        end
        if ~isempty(Ni_index)
            sum_Ni= sum_Dt(Ni_index);
            r_index=find(Ni_index==vacancy_atomID);
            sum_Ni(r_index)=[];
            MSD_Ni(i+1) =sum(sum_Ni)./length(Ni_index);
        end
        
    end
   
    t_set=tracks(:,1);
    MSD_t=tracks(:,5);
    disp(['dif= ',num2str(sum(MSD_t-MSD))])
    p0= polyfit(t_set, MSD_t,1);
    p1= polyfit(t_set, MSD_Ni,1);
    p2= polyfit(t_set, MSD_Co,1);
    p3= polyfit(t_set, MSD_Cr,1);
    x0=0:max(t_set)/1e3:max(t_set);
    y0= polyval(p0,x0);   
    y1= polyval(p1,x0);    
    y2= polyval(p2,x0);
    y3= polyval(p3,x0);
%     figure;
%     scatter(t_set,MSD_t,'.','MarkerEdgeColor','#1b9e77')
%     hold on
%     plot(x0,y0,'-','Color','#1b9e77','LineWidth',Lw);
%     scatter(t_set,MSD_Ni,'.','MarkerEdgeColor','#7570b3')
%     plot(x0,y1,'-','Color','#7570b3','LineWidth',Lw);
%     scatter(t_set,MSD_Co,'.','MarkerEdgeColor','#e6ab02')
%     plot(x0,y2,'-','Color','#e6ab02','LineWidth',Lw);
%     scatter(t_set,MSD_Cr,'.','MarkerEdgeColor','#d95f02')
%     plot(x0,y3,'-','Color','#d95f02','LineWidth',Lw);
%     hold off;
%     xlabel('time(s)')
%     ylabel('Å^2')
%     legend('NiCoCr','fit-NiCoCr','Ni','fit-Ni','Co','fit-Co','Cr','fit-Cr','box','off')
%     title([num2str(T),'K']);
%     set(gca,'FontName','Arial','FontSize',12,'FontWeight','bold','Linewidth',Lw)
    %       saveas(gcf,['Ni(',num2str((num_ratio-1)/10),')Fe(',num2str(1-(num_ratio-1)/10),').fig'])
    D(num_ratio)=p0(1)./6;
    D_Co(num_ratio)=p1(1)./6;
    D_Ni(num_ratio)=p2(1)./6;
    D_Cr(num_ratio)=p3(1)./6;
    a=norm(current_jump_vector);
    Ni_mig_num= length(find(mig_type_set==1));
    Co_mig_num= length(find(mig_type_set==2));
    Cr_mig_num= length(find(mig_type_set==3));
    f_tracer_all(num_ratio)= MSD_t(end)/((length(MSD_t)/N) * a^2);
    if  ~isempty(Co_index)
        f_tracer_Co(num_ratio)= MSD_Co(end)/( Co_mig_num/length(Co_index) * a^2);
    end
    if  ~isempty(Ni_index)
        f_tracer_Ni(num_ratio)= MSD_Ni(end)/(Ni_mig_num/length(Ni_index) * a^2);
    end
    if  ~isempty(Cr_index)
        f_tracer_Cr(num_ratio)= MSD_Cr(end)/( Cr_mig_num/length(Cr_index) * a^2);
    end
    %     save(['ratio',num2str(num_ratio),'.mat'])
    per_coords_current=per_coords_set{num_ratio};
    dif_m=abs(R_r_vacancy- per_coords_current);
    dif_m=round(dif_m.*1e5)./1e5;
    mod_dif_m= mod(dif_m,Len);
    disp('check the per_coords(remove 1 frame):')
    disp(sum(sum(mod_dif_m)))  %%
    
        
end
save(['restart_SRO_T_',num2str(T),'K',num2str(nn),'nn.mat'],...
    'SRO_cell', 'strcture_cell', 'lattice_constant_set', 'ratio', 'T');