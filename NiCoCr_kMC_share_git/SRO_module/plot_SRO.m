function plot_SRO(x_axis, SRO_set,legend_set,xlabel_set,ylabel_set,title_set)
color_set={'#0072bd';
    '#d95319';
    '#edb120';
    '#7e2f8e';
    '#77ac30';
    '#4dbeee';
    '#a6cee3';
    '#1f78b4';
    '#b2df8a';
    '#33a02c';
    '#fb9a99';
    '#e31a1c';
    '#fdbf6f';
    '#ff7f00';
    '#cab2d6'};
Marker_set={'o','s','*','d','p','v','^','<','>','+'};

% SRO_count_set=[zeros(9,1),SRO];
figure;
MS=8;
LW=2;
x=x_axis;
y=SRO_set;
% x=[0:1:10].*10;
% y=SRO_count_set;
plot(x,y(1,:)./2,'MarkerSize',MS, 'Color',color_set{1},'Linewidth',LW,'Marker',Marker_set{1});
hold on;
plot(x,y(6,:),'MarkerSize',MS, 'Color',color_set{2},'Linewidth',LW,'Marker',Marker_set{6});
plot(x,y(2,:),'MarkerSize',MS, 'Color',color_set{3},'Linewidth',LW,'Marker',Marker_set{2});
plot(x,y(3,:),'MarkerSize',MS, 'Color',color_set{4},'Linewidth',LW,'Marker',Marker_set{10});
% plot(x,y(3,:),'MarkerSize',MS, 'Color',color_set{4},'Linewidth',LW,'Marker',Marker_set{4});
% plot(x,y(4,:),'MarkerSize',MS, 'Color',color_set{5},'Linewidth',LW);
plot(x,y(5,:)./2,'MarkerSize',MS, 'Color',color_set{5},'Linewidth',LW,'Marker',Marker_set{7});

% plot(x,y(7,:),'MarkerSize',MS, 'Color',color_set{8},'Linewidth',LW);
% plot(x,y(8,:),'MarkerSize',MS, 'Color',color_set{9},'Linewidth',LW);
plot(x,y(9,:)./2,'MarkerSize',MS, 'Color',color_set{6},'Linewidth',LW,'Marker',Marker_set{4});
% plot(plot_set,Sro_set(plot_set,1),'Marker',Marker_set{1},'MarkerSize',MS, 'Color',color_set{1},'Linewidth',LW,'Linestyle',':');
% hold on;plot(x,y(8,:),'MarkerSize',MS, 'Color',color_set{9},'Linewidth',LW);
% plot(plot_set,Sro_set(plot_set,2),'Marker',Marker_set{2},'MarkerSize',MS, 'Color',color_set{2},'Linewidth',LW,'Linestyle','--');
% plot(plot_set,Sro_set(plot_set,3),'Marker',Marker_set{3},'MarkerSize',MS, 'Color',color_set{3},'Linewidth',LW,'Linestyle',':');
% plot(plot_set,Sro_set(plot_set,4),'Marker',Marker_set{4},'MarkerSize',MS, 'Color',color_set{4},'Linewidth',LW,'Linestyle','--');
% plot(plot_set,Sro_set(plot_set,5),'Marker',Marker_set{5},'MarkerSize',MS, 'Color',color_set{5},'Linewidth',LW,'Linestyle',':');
% plot(plot_set,Sro_set(plot_set,6),'Marker',Marker_set{6},'MarkerSize',MS, 'Color',color_set{6},'Linewidth',LW,'Linestyle','--');
% plot(plot_set,Sro_set(plot_set,7),'Marker',Marker_set{7},'MarkerSize',MS, 'Color',color_set{7},'Linewidth',LW,'Linestyle',':');
% plot(plot_set,Sro_set(plot_set,8),'Marker',Marker_set{8},'MarkerSize',MS, 'Color',color_set{8},'Linewidth',LW,'Linestyle','--');
% plot(plot_set,Sro_set(plot_set,9),'Marker',Marker_set{9},'MarkerSize',MS, 'Color',color_set{9},'Linewidth',LW,'Linestyle',':');
% legend('Ni-Ni','Ni-Co','Ni-Cr','Co-Co','Co-Cr','Cr-Cr','box','off')
% xlabel('Time (ns)')
% ylabel('WC-SRO')
% title('1350K')
legend(legend_set,'box','off')
xlabel(xlabel_set)
ylabel(ylabel_set)
title(title_set)
ylim([-0.8,0.8]);
grid on;
set(gca, 'FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',2)
% saveas(gcf,'SRO_change.fig')