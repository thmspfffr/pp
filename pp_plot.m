%%

clear
restoredefaultpath

addpath ~/pconn/matlab/
load /home/gnolte/meth/templates/mri.mat
load /home/gnolte/meth/templates/sa_template.mat
addpath ~/Documents/MATLAB/fieldtrip-20190224/
ft_defaults



addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
v = 1;

[plt_gla,plt_hh,plt_mue,plt_all]=pp_load_results(v);

colors = cbrewer('qual', 'Set3', 10,'pchip'); 
colors = colors(4:6,:);

%% PLOT MUTUAL INFORMATION
freqoi=2.^(1:(1/4):7); 

figure_w

subplot(2,3,1);
imagesc(nanmean(plt_gla.sens_mi_ord,3),[0 .02])
colormap(plasma); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,2);


subplot(2,3,3);


%% PLOT CORRELATION IN SENSOR SPACE - SORTED FROM ANTERIOR TO POSTERIOR
freqoi=2.^(1:(1/4):7); 

figure_w

subplot(2,3,1);
imagesc(nanmean(plt_gla.corr_sens_ord(3:end,:,:),3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,2);
imagesc(nanmean(plt_hh.corr_sens_ord(3:end,:,:),3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,3);
imagesc(nanmean(plt_mue.corr_sens_ord(3:end,:,:),3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

pooled = cat(3,plt_hh.corr_sens_ord,plt_gla.corr_sens_ord,plt_mue.corr_sens_ord);

subplot(2,3,4);
imagesc(nanmean(pooled(3:end,:,:),3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3);
h = p<fdr1(p(:),0.05,0);

subplot(2,3,5);
imagesc(nanmean(pooled(3:end,:,:),3).*h(3:end,:),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_anterior_post_v%d.pdf',v))

%% LOW AND HIGH FREQS


lf = nanmean(squeeze(mean(plt_gla.corr_sens_ord(3:end,1:5,:),2)),2);
hf = nanmean(squeeze(mean(plt_gla.corr_sens_ord(3:end,21:25,:),2)),2);

[r_gla, p_gla] = corr(nanmean(lf,2),nanmean(hf,2));
    

figure_w;

subplot(2,3,1)
scatter(zscore(nanmean(lf,2)),zscore(nanmean(hf,2)),30,'markeredgecolor','w','markerfacecolor',colors(1,:)); axis square
lsline; xlabel('Correlation (2-4 Hz)'); ylabel('Correlation (64-128 Hz)')
axis([-3 3 -3 3]); tp_editplots
set(gca,'xtick',[-3 0 3],'xticklabel',[-3 0 3])
set(gca,'ytick',[-3 0 3],'yticklabel',[-3 0 3])
text(-1.5,-2.5,sprintf('r = %.2f | p < %.2f',r_gla,p_gla),'fontsize',7)

lf = nanmean(squeeze(mean(plt_hh.corr_sens_ord(3:end,1:5,:),2)),2);
hf = nanmean(squeeze(mean(plt_hh.corr_sens_ord(3:end,21:25,:),2)),2);

[r_hh, p_hh] = corr(nanmean(lf,2),nanmean(hf,2));

subplot(2,3,2)
scatter(zscore(nanmean(lf,2)),zscore(nanmean(hf,2)),30,'markeredgecolor','w','markerfacecolor',colors(2,:)); axis square
lsline; xlabel('Correlation (2-4 Hz)')
axis([-3 3 -3 3]); tp_editplots
set(gca,'xtick',[-3 0 3],'xticklabel',[-3 0 3])
set(gca,'ytick',[-3 0 3],'yticklabel',[-3 0 3])
text(-1.5,-2.5,sprintf('r = %.2f | p < %.2f',r_hh,p_hh),'fontsize',7)

lf = nanmean(squeeze(mean(plt_mue.corr_sens_ord(3:end,1:5,:),2)),2);
hf = nanmean(squeeze(mean(plt_mue.corr_sens_ord(3:end,21:25,:),2)),2);

[r_mue, p_mue] = corr(nanmean(lf,2),nanmean(hf,2));

subplot(2,3,3)
scatter(zscore(nanmean(lf,2)),zscore(nanmean(hf,2)),30,'markeredgecolor','w','markerfacecolor',colors(3,:)); axis square
lsline; xlabel('Correlation (2-4 Hz)')
axis([-3 3 -3 3]); tp_editplots
set(gca,'xtick',[-3 0 3],'xticklabel',[-3 0 3])
set(gca,'ytick',[-3 0 3],'yticklabel',[-3 0 3])
text(-1.5,-2.5,sprintf('r = %.2f | p < %.2f',r_mue,p_mue),'fontsize',7)

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_correlation_lowhigh_v%d.pdf',v))


%% PLOT AVERAGES ACROSS SPACE, HAMBURG, GLASGOW, ALL

% -----------------
% SENSOR SPACE - POWER
% -----------------
freqoi=2.^(1:(1/4):7); 

tot_size = size(plt_mue.corr_sens,3)+size(plt_gla.corr_sens,3)+size(plt_hh.corr_sens,3);

par_gla = mean(log10(nanmean(plt_gla.pow_sens,1)),3);
par_hh  = mean(log10(nanmean(plt_hh.pow_sens,1)),3);
par_mue = nanmean(log10(nanmean(plt_mue.pow_sens,1)),3);

std_gla = std(log10(nanmean(plt_gla.pow_sens,1)),[],3)/sqrt(size(plt_gla.pow_sens,3));
std_hh  = std(log10(nanmean(plt_hh.pow_sens,1)),[],3)/sqrt(size(plt_hh.pow_sens,3));
std_mue = nanstd(log10(nanmean(plt_mue.pow_sens,1)),[],3)/sqrt(size(plt_mue.pow_sens,3));

par_all = mean(cat(2,log10(squeeze(nanmean(plt_gla.pow_sens,1))),log10(squeeze(nanmean(plt_hh.pow_sens,1))),log10(squeeze(nanmean(plt_mue.pow_sens,1)))),2);
std_all = std(cat(2,log10(squeeze(nanmean(plt_gla.pow_sens,1))),log10(squeeze(nanmean(plt_hh.pow_sens,1))),log10(squeeze(nanmean(plt_mue.pow_sens,1)))),[],2)/sqrt(tot_size);

figure_w
subplot(4,3,1)
shadedErrorBar(log10(freqoi),par_gla,std_gla,{'color',colors(1,:)})
axis([.3 2.11 -24 -21])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,2)
shadedErrorBar(log10(freqoi),par_hh,std_hh,{'color',colors(2,:)})
axis([.3 2.11 -24 -21])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,3)
shadedErrorBar(log10(freqoi),par_mue,std_mue,{'color',colors(3,:)})
axis([.3 2.11 -24 -21])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,4)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -24 -21])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_pow_lineplots_v%d.pdf',v))


figure_w

% -----------------
% SENSOR SPACE - CORRELATIONS
% -----------------
par_gla = mean(nanmean(plt_gla.corr_sens,1),3);
par_hh  = mean(nanmean(plt_hh.corr_sens,1),3);
par_mue = mean(nanmean(plt_mue.corr_sens,1),3);

std_gla = std(nanmean(plt_gla.corr_sens,1),[],3)/sqrt(size(plt_gla.corr_sens,3));
std_hh  = std(nanmean(plt_hh.corr_sens,1),[],3)/sqrt(size(plt_hh.corr_sens,3));
std_mue = std(nanmean(plt_mue.corr_sens,1),[],3)/sqrt(size(plt_mue.corr_sens,3));

par_all = mean(cat(2,squeeze(nanmean(plt_gla.corr_sens,1)),squeeze(nanmean(plt_hh.corr_sens)),squeeze(nanmean(plt_mue.corr_sens))),2);
std_all = std(cat(2,squeeze(nanmean(plt_gla.corr_sens,1)),squeeze(nanmean(plt_hh.corr_sens)),squeeze(nanmean(plt_mue.corr_sens))),[],2)/sqrt(50);

subplot(4,3,1)
shadedErrorBar(log10(freqoi),par_gla,std_gla,{'color',colors(1,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,2)
shadedErrorBar(log10(freqoi),par_hh,std_hh,{'color',colors(2,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,3)
shadedErrorBar(log10(freqoi),par_mue,std_mue,{'color',colors(3,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,4)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_corr_lineplots_v%d.pdf',v))

% -----------------
% SOURCE SPACE: NORMAL PUPIL 
% -----------------

figure_w

par_gla = mean(nanmean(plt_gla.corr_src,1),3);
par_hh  = mean(nanmean(plt_hh.corr_src,1),3);
par_mue = mean(nanmean(plt_mue.corr_src,1),3);

std_gla = std(nanmean(plt_gla.corr_src,1),[],3)/sqrt(size(plt_gla.corr_src,3));
std_hh  = std(nanmean(plt_hh.corr_src,1),[],3)/sqrt(size(plt_hh.corr_src,3));
std_mue = std(nanmean(plt_mue.corr_src,1),[],3)/sqrt(size(plt_mue.corr_src,3));

par_all = nanmean(mean(cat(3,plt_gla.corr_src,plt_hh.corr_src,plt_mue.corr_src),1),3);
std_all = std(nanmean(cat(3,plt_gla.corr_src,plt_hh.corr_src,plt_mue.corr_src),1),[],3)/sqrt(tot_size);

subplot(4,3,1)
shadedErrorBar(log10(freqoi),par_gla,std_gla,{'color',colors(1,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,2)
shadedErrorBar(log10(freqoi),par_hh,std_hh,{'color',colors(2,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,3)
shadedErrorBar(log10(freqoi),par_mue,std_mue,{'color',colors(3,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,4)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_lineplots_v%d.pdf',v))

% --------------------------------------------
% SOURCE SPACE AVERAGES: PUPIL DERIVATIVE
% --------------------------------------------

par_gla = mean(nanmean(plt_gla.corr_src_df,1),3);
par_hh = mean(mean(plt_hh.corr_src_df,1),3);
par_mue = mean(nanmean(plt_mue.corr_src_df,1),3);

std_gla = std(nanmean(plt_gla.corr_src_df,1),[],3)/sqrt(size(plt_gla.corr_src_df,3));
std_hh  = std(nanmean(plt_hh.corr_src_df,1),[],3)/sqrt(size(plt_hh.corr_src_df,2));
std_mue  = std(nanmean(plt_mue.corr_src_df,1),[],3)/sqrt(size(plt_mue.corr_src_df,2));

par_all = nanmean(mean(cat(3,plt_gla.corr_src_df,plt_hh.corr_src_df,plt_mue.corr_src_df),1),3);
std_all = std(nanmean(cat(3,plt_gla.corr_src_df,plt_hh.corr_src_df,plt_mue.corr_src_df),1),[],3)/sqrt(tot_size);

figure_w;
subplot(4,3,1)
shadedErrorBar(log10(freqoi),par_gla,std_gla,{'color',colors(1,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
% title('Derivative: GLA')

subplot(4,3,2)
shadedErrorBar(log10(freqoi),par_hh,std_hh,{'color',colors(2,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
% title('Derivative: HH')

subplot(4,3,3)
shadedErrorBar(log10(freqoi),par_mue,std_mue,{'color',colors(3,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
% title('Derivative: MUE')

subplot(4,3,4)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
% title('Derivative: POOLED')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_derivative_lineplots_v%d.pdf',v))

%% PLOT SOURCE SPACE 
figure_w
subplot(2,4,1)
[h,p]=ttest(plt_gla.corr_src_BNA,zeros(size(plt_gla.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_gla.corr_src_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,2)
[h,p]=ttest(plt_hh.corr_src_BNA,zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_hh.corr_src_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,3)
[h,p]=ttest(plt_mue.corr_src_BNA,zeros(size(plt_mue.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_mue.corr_src_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,4)

pooled= cat(3,plt_hh.corr_src_BNA,plt_gla.corr_src_BNA,plt_mue.corr_src_BNA); 

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,0);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,5)
[h,p]=ttest(plt_gla.corr_src_df_BNA,zeros(size(plt_gla.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_gla.corr_src_df_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,6)
[h,p]=ttest(plt_hh.corr_src_df_BNA,zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_hh.corr_src_df_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,7)
[h,p]=ttest(plt_mue.corr_src_df_BNA,zeros(size(plt_mue.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_mue.corr_src_df_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,8)

pooled= cat(3,plt_hh.corr_src_df_BNA,plt_gla.corr_src_df_BNA,plt_mue.corr_src_df_BNA); 

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,0);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_correlations_BNA_v%d.pdf',v))

%% PLOT REGIONS OF INTEREST, SENSORY AREAS
load ~/pp/proc/pp_atlas_BNA.mat
M1 = [-42 -26 54; 38 -32 48];
V1 = [-20 -86 18; 16 -80 26];
A1 = [-54 -22 10; 52 -24 12];

for i= 1 : size(BNA.grid_5mm,1)
  
  d_A1_l(i) = sqrt((BNA.grid_5mm(i,1)-A1(1,1))^2 + (BNA.grid_5mm(i,2)-A1(1,2))^2 + (BNA.grid_5mm(i,3)-A1(1,3))^2);
  d_A1_r(i) = sqrt((BNA.grid_5mm(i,1)-A1(2,1))^2 + (BNA.grid_5mm(i,2)-A1(2,2))^2 + (BNA.grid_5mm(i,3)-A1(2,3))^2);
  d_V1_l(i) = sqrt((BNA.grid_5mm(i,1)-V1(1,1))^2 + (BNA.grid_5mm(i,2)-V1(1,2))^2 + (BNA.grid_5mm(i,3)-V1(1,3))^2);
  d_V1_r(i) = sqrt((BNA.grid_5mm(i,1)-V1(2,1))^2 + (BNA.grid_5mm(i,2)-V1(2,2))^2 + (BNA.grid_5mm(i,3)-V1(2,3))^2);
  d_M1_l(i) = sqrt((BNA.grid_5mm(i,1)-M1(1,1))^2 + (BNA.grid_5mm(i,2)-M1(1,2))^2 + (BNA.grid_5mm(i,3)-M1(1,3))^2);
  d_M1_r(i) = sqrt((BNA.grid_5mm(i,1)-M1(2,1))^2 + (BNA.grid_5mm(i,2)-M1(2,2))^2 + (BNA.grid_5mm(i,3)-M1(2,3))^2);
  
end
[~,A1_l]=min(d_A1_l);
[~,A1_r]=min(d_A1_r);
[~,V1_l]=min(d_V1_l);
[~,V1_r]=min(d_V1_r);
[~,M1_l]=min(d_M1_l);
[~,M1_r]=min(d_M1_r);

figure_w;
subplot(4,3,1); hold on
r_A1_lr = corr(nanmean(plt_all.corr_src(A1_l,:,:),3)',nanmean(plt_all.corr_src(A1_r,:,:),3)');
sig_A1_lr = mean(plt_all.corr_src([A1_l A1_r],:,:),1);
[~,p]=ttest(sig_A1_lr,zeros(size(sig_A1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('A1: r = %.3f',r_A1_lr))

subplot(4,3,2); hold on
r_V1_lr = corr(nanmean(plt_all.corr_src(V1_l,:,:),3)',nanmean(plt_all.corr_src(V1_r,:,:),3)');
sig_V1_lr = mean(plt_all.corr_src([V1_l V1_r],:,:),1);
[~,p]=ttest(sig_V1_lr,zeros(size(sig_V1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,3); hold on
r_M1_lr = corr(nanmean(plt_all.corr_src(M1_l,:,:),3)',nanmean(plt_all.corr_src(M1_r,:,:),3)');
sig_M1_lr = mean(plt_all.corr_src([M1_l M1_r],:,:),1);
[~,p]=ttest(sig_M1_lr,zeros(size(sig_M1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('SOM: r = %.3f',r_M1_lr))

subplot(4,3,4); hold on
r_A1_lr = corr(nanmean(plt_all.corr_src_df(A1_l,:,:),3)',nanmean(plt_all.corr_src_df(A1_r,:,:),3)');
sig_A1_lr = mean(plt_all.corr_src_df([A1_l A1_r],:,:),1);
[~,p]=ttest(sig_A1_lr,zeros(size(sig_A1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('A1: r = %.3f',r_A1_lr))

subplot(4,3,5); hold on
r_V1_lr = corr(nanmean(plt_all.corr_src_df(V1_l,:,:),3)',nanmean(plt_all.corr_src_df(V1_r,:,:),3)');
sig_V1_lr = mean(plt_all.corr_src_df([V1_l V1_r],:,:),1);
[~,p]=ttest(sig_V1_lr,zeros(size(sig_V1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,6); hold on
r_M1_lr = corr(nanmean(plt_all.corr_src_df(M1_l,:,:),3)',nanmean(plt_all.corr_src_df(M1_r,:,:),3)');
sig_M1_lr = mean(plt_all.corr_src_df([M1_l M1_r],:,:),1);
[~,p]=ttest(sig_M1_lr,zeros(size(sig_M1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('SOM: r = %.3f',r_M1_lr))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_v%d.pdf',v))

%% PLOT SOURCE MAPS: PUPIL
% load /home/gnolte/meth/templates/mri.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

for ifoi = 1:25
% ifoi = 14;

[h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
h=p<(fdr1(p(:),0.05,0));
par=nanmean(plt_all.corr_src(:,ifoi,:),3).*h;
% project onto fine grid
par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);

clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.95;
para.colorbar= 0;

tp_showmri_transp(mri,para,[sa_template.grid_fine par]); f=get(gcf)

move_x = [0:0.08:0.32];
move_y = [0.08:-0.02:0];

for i=4:-1:1
    for j = 0:3
        if j ~= 4
            f.Children(i*4-j).Position(1)=f.Children(i*4-j).Position(1)-move_x(j+1);
        end
        f.Children(i*4-j).Position(2)=f.Children(i*4-j).Position(2)+move_y(i);
    end
end


set(gcf,'renderer','painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_avg_f%d_v%d.tiff',ifoi,v))

end

%% PLOT SOURCE MAPS: PUPIL DERIVATIVE
% load /home/gnolte/meth/templates/mri.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

for ifoi = 1:25
% ifoi = 14;

[h,p] = ttest(plt_all.corr_src_df(:,ifoi,:),zeros(size(plt_all.corr_src_df(:,ifoi,:))),'dim',3);
h=p<(fdr1(p(:),0.05,0));
par=nanmean(plt_all.corr_src_df(:,ifoi,:),3).*h;
% project onto fine grid
par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);

clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.95;
para.colorbar= 0;


tp_showmri_transp(mri,para,[sa_template.grid_fine par]); f=get(gcf)

move_x = [0:0.08:0.32];
move_y = [0.08:-0.02:0];

for i=4:-1:1
    for j = 0:3
        if j ~= 4
            f.Children(i*4-j).Position(1)=f.Children(i*4-j).Position(1)-move_x(j+1);
        end
        f.Children(i*4-j).Position(2)=f.Children(i*4-j).Position(2)+move_y(i);
    end
end

set(gcf,'renderer','painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_df_sourcemap_avg_f%d_v%d.tiff',ifoi,v))

end

error('!')
%% PLOT COUNTING/TASK RESULTS 

figure_w

subplot(2,3,1)
imagesc(nanmean(plt_hh.corr_src_BNA_cnt,3),[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,2)
imagesc(nanmean(plt_hh.corr_src_BNA_cnt,3)-nanmean(plt_hh.corr_src_BNA,3),[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,3)
[h,p]=ttest(plt_hh.corr_src_BNA_cnt,zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc((nanmean(plt_hh.corr_src_BNA_cnt,3)-nanmean(plt_hh.corr_src_BNA,3)).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,4)
[h,p]=ttest(plt_hh.corr_src_BNA_cnt,plt_hh.corr_src_BNA,'dim',3); 
h=p<0.05;
imagesc((nanmean(plt_hh.corr_src_BNA_cnt,3)-nanmean(plt_hh.corr_src_BNA,3)).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_correlations_BNA_counting_v%d.pdf',v))
%% PLOT COUNTING: REGIONS OF INTEREST, SENSORY AREAS

load ~/pp/proc/pp_atlas_BNA.mat
M1 = [-42 -26 54; 38 -32 48];
V1 = [-20 -86 18; 16 -80 26];
A1 = [-54 -22 10; 52 -24 12];

for i= 1 : size(BNA.grid_5mm,1)
  
  d_A1_l(i) = sqrt((BNA.grid_5mm(i,1)-A1(1,1))^2 + (BNA.grid_5mm(i,2)-A1(1,2))^2 + (BNA.grid_5mm(i,3)-A1(1,3))^2);
  d_A1_r(i) = sqrt((BNA.grid_5mm(i,1)-A1(2,1))^2 + (BNA.grid_5mm(i,2)-A1(2,2))^2 + (BNA.grid_5mm(i,3)-A1(2,3))^2);
  d_V1_l(i) = sqrt((BNA.grid_5mm(i,1)-V1(1,1))^2 + (BNA.grid_5mm(i,2)-V1(1,2))^2 + (BNA.grid_5mm(i,3)-V1(1,3))^2);
  d_V1_r(i) = sqrt((BNA.grid_5mm(i,1)-V1(2,1))^2 + (BNA.grid_5mm(i,2)-V1(2,2))^2 + (BNA.grid_5mm(i,3)-V1(2,3))^2);
  d_M1_l(i) = sqrt((BNA.grid_5mm(i,1)-M1(1,1))^2 + (BNA.grid_5mm(i,2)-M1(1,2))^2 + (BNA.grid_5mm(i,3)-M1(1,3))^2);
  d_M1_r(i) = sqrt((BNA.grid_5mm(i,1)-M1(2,1))^2 + (BNA.grid_5mm(i,2)-M1(2,2))^2 + (BNA.grid_5mm(i,3)-M1(2,3))^2);
  
end
[~,A1_l]=min(d_A1_l);
[~,A1_r]=min(d_A1_r);
[~,V1_l]=min(d_V1_l);
[~,V1_r]=min(d_V1_r);
[~,M1_l]=min(d_M1_l);
[~,M1_r]=min(d_M1_r);

figure_w;
subplot(4,3,1); hold on
r_A1_lr = corr(nanmean(plt_hh.corr_src_cnt(A1_l,:,:),3)',nanmean(plt_hh.corr_src_cnt(A1_r,:,:),3)');
sig_A1_lr = mean(plt_hh.corr_src_cnt([A1_l A1_r],:,:),1);
[~,p]=ttest(sig_A1_lr,zeros(size(sig_A1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('A1: r = %.3f',r_A1_lr))

subplot(4,3,2); hold on
r_V1_lr = corr(nanmean(plt_hh.corr_src_cnt(V1_l,:,:),3)',nanmean(plt_hh.corr_src_cnt(V1_r,:,:),3)');
sig_V1_lr = mean(plt_hh.corr_src_cnt([V1_l V1_r],:,:),1);
[~,p]=ttest(sig_V1_lr,zeros(size(sig_V1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,3); hold on
r_M1_lr = corr(nanmean(plt_hh.corr_src_cnt(M1_l,:,:),3)',nanmean(plt_hh.corr_src_cnt(M1_r,:,:),3)');
sig_M1_lr = mean(plt_hh.corr_src_cnt([M1_l M1_r],:,:),1);
[~,p]=ttest(sig_M1_lr,zeros(size(sig_M1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
plot(log10(freqoi(h)),nanmean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); title(sprintf('SOM: r = %.3f',r_M1_lr))


subplot(4,3,4); hold on
sig_A1_lr = mean(plt_hh.corr_src_cnt([A1_l A1_r],:,:),1)-mean(plt_hh.corr_src([A1_l A1_r],:,:),1);
[~,p]=ttest(sig_A1_lr,zeros(size(sig_A1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
plot(log10(freqoi(h)),mean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); 

subplot(4,3,5); hold on
sig_V1_lr = mean(plt_hh.corr_src_cnt([V1_l V1_r],:,:),1)-mean(plt_hh.corr_src([V1_l V1_r],:,:),1);
[~,p]=ttest(sig_V1_lr,zeros(size(sig_V1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
plot(log10(freqoi(h)),mean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); 

subplot(4,3,6); hold on
sig_M1_lr = mean(plt_hh.corr_src_cnt([M1_l M1_r],:,:),1)-mean(plt_hh.corr_src([M1_l M1_r],:,:),1);
[~,p]=ttest(sig_M1_lr,zeros(size(sig_M1_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
plot(log10(freqoi(h)),mean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); 


print(gcf,'-dpdf',sprintf('~/pp/plots/pp_cnt_src_ROIs_v%d.pdf',v))

%% PLOT SOURCE

addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% ifoi = 5
ifoi = 10:12;

[h,p] = ttest(nanmean(plt_hh.corr_src_cnt(:,ifoi,:),2),nanmean(plt_hh.corr_src(:,ifoi,:),2),'dim',3);

% h=p<(fdr1(p(:),0.2,0));
par=(nanmean(nanmean(plt_hh.corr_src_cnt(:,ifoi,:),2),3)-nanmean(nanmean(plt_hh.corr_src(:,ifoi,:),2),3)).*h;

clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.75;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par])
set(gcf,'renderer','painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_cnt_src_corr_sourcemap_avg_f%df%d_v%d.tiff',min(ifoi),max(ifoi),v))

  
%% PLOT CROSS FREQUENCY
meg_f1 = [2 3 4 6]; pup_f1 = [9 10 11 12 13 14];
meg_f2 = [9 10 11]; pup_f2 = [6 7 8 9 10];
meg_f3 = [13 14 15]; pup_f3 = [6 7 8 9 10 11];
meg_f4 = [10:12];  pup_f4 = [1:5];
meg_f5 = [21:25];  pup_f5 = [2 3 4 5 6];
meg_f6 = [9:13];  pup_f6 = [12:17];

pooled = cat(4,plt_hh.cf_corr,plt_gla.cf_corr);
pupil_freqoi(:,1) = 2.^(-9:(0.5):1);
pupil_freqoi(:,2) = 2.^(-8:(0.5):2);


figure_w;

subplot(2,3,1)
imagesc(squeeze(nanmean(nanmean(plt_gla.cf_corr,4),1))',[-0.02 0.02]);
set(gca,'ydir','normal'); axis square
set(gca,'tickdir','out','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
set(gca,'tickdir','out','ytick',1:2:21,'yticklabel',num2str((round((mean(pupil_freqoi(1:2:21,:),2)*10000))/10000)))
tp_editplots
xlabel('MEG frequency [Hz]'); ylabel('Pupil frequency [Hz]')
colormap(cmap)

subplot(2,3,2)
imagesc(squeeze(nanmean(nanmean(plt_hh.cf_corr,4),1))',[-0.02 0.02]);
set(gca,'ydir','normal'); axis square
set(gca,'tickdir','out','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
set(gca,'tickdir','out','ytick',1:2:21,'yticklabel',num2str((round((mean(pupil_freqoi(1:2:21,:),2)*10000))/10000)))
tp_editplots
xlabel('MEG frequency [Hz]'); ylabel('Pupil frequency [Hz]')
colormap(cmap)

subplot(2,3,3)
imagesc(squeeze(nanmean(nanmean(pooled,4),1))',[-0.02 0.02]);
set(gca,'ydir','normal'); axis square
set(gca,'tickdir','out','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
set(gca,'tickdir','out','ytick',1:2:21,'yticklabel',num2str((round((mean(pupil_freqoi(1:2:21,:),2)*10000))/10000)))

tp_drawsquare(min(meg_f1), max(meg_f1), min(pup_f1), max(pup_f1))
tp_drawsquare(min(meg_f2), max(meg_f2), min(pup_f2), max(pup_f2))
tp_drawsquare(min(meg_f3), max(meg_f3), min(pup_f3), max(pup_f3))
tp_drawsquare(min(meg_f4), max(meg_f4), min(pup_f4), max(pup_f4))
tp_drawsquare(min(meg_f5), max(meg_f5), min(pup_f5), max(pup_f5))
tp_drawsquare(min(meg_f6), max(meg_f6), min(pup_f6), max(pup_f6))

tp_editplots
xlabel('MEG frequency [Hz]'); ylabel('Pupil frequency [Hz]')
colormap(cmap)

subplot(2,3,4)
[h,p]=ttest(squeeze(nanmean(pooled,1)),plt_hhzeros(size(squeeze(nanmean(pooled,1)))),'dim',3);
h=p<0.05;
imagesc(squeeze(nanmean(nanmean(pooled,4),1))'.*h',[-0.02 0.02]);
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
set(gca,'ytick',1:2:21,'yticklabel','')
tp_editplots; axis square
xlabel('MEG frequency [Hz]'); 
colormap(cmap)

v=3
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_crossfreq_v%d.pdf',v))

%% IMAGE AGAIN



addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

ifreq = 6
% [h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
% h=p<(fdr1(p(:),0.05,0));
par=nanmean(nanmean(nanmean(plt_hh.cf_corr(:,eval(sprintf('meg_f%d',ifreq)),eval(sprintf('pup_f%d',ifreq)),:),4),3),2);

par_stats=squeeze(nanmean(nanmean(plt_hh.cf_corr(:,eval(sprintf('meg_f%d',ifreq)),eval(sprintf('pup_f%d',ifreq)),:),3),2));
% [h,p]=ttest(par_stats,zeros(size(par_stats)),'dim',2); h = p<(fdr1(p(:),0.05,0));
par=par;

clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.75;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par])
set(gcf,'renderer','painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_crossfreq_f%d_v%d.tiff',ifreq,v))



%% FRACTION OF POSITIVE/NEGATIVE CORRELATIONS
h=ttest(pooled,zeros(size(pooled)),'dim',4);
par = nanmean(pooled,4);
par_pos = 100*(squeeze(sum((h>0) & (par>0)))./8799);
par_neg = 100*(squeeze(sum((h>0) & (par<0)))./8799);
% permtest
% idx=randi(2,[50 1000])-1;

% all_dat = cat(4,pooled,zeros(size(pooled)));
% nsubj = size(pooled,4);
% for iperm = 1 : 1000
%   iperm
%   clear par
%   par = all_dat(:,:,:,(1:50)'+nsubj.*idx(:,iperm));
% 
%   h=ttest(par,zeros(size(pooled)),'dim',4);
%   par_perm = nanmean(par,4);
%   par_pos_perm(:,:,iperm) = 100*(squeeze(sum((h>0) & (par_perm>0)))./8799);
%   par_neg_perm(:,:,iperm) = 100*(squeeze(sum((h>0) & (par_perm<0)))./8799);
% 
% end

mask_pos = (1-sum(par_pos>par_pos_perm,3)/1000)<0.05;
mask_neg = (1-sum(par_neg>par_neg_perm,3)/1000)<0.05;

par = zeros(size(par_pos));
par(mask_pos)=par_pos(mask_pos);
par(mask_neg)=-par_neg(mask_neg);


k=subplot(2,3,4)
imagesc(par',[-75 75]);
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
set(gca,'ytick',1:2:21,'yticklabel','')
tp_editplots; axis square
xlabel('MEG frequency [Hz]'); 
colormap(cmap)

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_crossfreq_v%d.pdf',v))

%% RECEPTORS: RAW PUPIL (i.e., not derivative)
load ~/receptormaps/gene_values.mat
orig = [46 64 37];
locs_mni=2*(locs-repmat(orig,[size(locs,1) 1]))/10;
isabs=1;
  
dd = 1;
for irecep = 1 : 54

  freqoi=2.^(1:(1/4):7);
  par = nan(8799,1);
  par(ismember(BNA.tissue_5mm,[1:2:246]))=spatfiltergauss(nanmean(gdat(:,irecep),2),locs_mni,dd,BNA.grid_5mm(ismember(BNA.tissue_5mm,[1:2:246]),:)/10);
  clear r p par2 par1
  for isubj = 1 : 50
    isubj
    k=0;
    pp = plt_all.corr_src(:,:,isubj);
    pp_df = plt_all.corr_src_df(:,:,isubj);
    for igrid = 1:2:246
      k=k+1;
      par1(k) = mean(par(BNA.tissue_5mm == igrid));
    end
    k = 0;
    for igrid = 1 : 246
      k= k+1;
      tmp(k,:) = mean(pp(BNA.tissue_5mm == igrid,:));
      tmp_df(k,:) = mean(pp_df(BNA.tissue_5mm == igrid,:));
    end
    par2=(tmp(1:2:246,:)+tmp(2:2:246,:))/2;
    par2_df=(tmp_df(1:2:246,:)+tmp_df(2:2:246,:))/2;
    %     par2(par2<0)=0;
    if isabs==1
      [r(:,isubj),p(:,isubj)] = corr(par1(:),abs(par2),'type','pearson');
      [r_df(:,isubj),p_df(:,isubj)] = corr(par1(:),abs(par2_df),'type','pearson');
    else
      [r(:,isubj),p(:,isubj)] = corr(par1(:),par2,'type','pearson');
      [r_df(:,isubj),p_df(:,isubj)] = corr(par1(:),par2_df,'type','pearson');
    end
  end
  

  [~,p]=ttest(r,zeros(size(r)),'dim',2);
  [~,p_df]=ttest(r_df,zeros(size(r_df)),'dim',2);

  if irecep == 1
    figure_w; 
    subplot(6,5,irecep); hold on
  elseif irecep ~= 1 && irecep<28
    subplot(6,5,irecep); hold on
  end
  
  if irecep == 28
    set(gcf,'Position',[50 50 800 1200])
    print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_receptors_abs%d_dd%s_nr1_v%d.pdf',isabs,num2str(dd),v))
    close
    figure_w; 
    subplot(6,5,irecep-27); hold on
  elseif irecep ~= 28 && irecep>27
    subplot(6,5,irecep-27); hold on
  end
    
  
  title(sprintf('%s',lab{irecep}))

  shadedErrorBar(1:25,mean(r,2),std(r,[],2)/sqrt(50),'k'); hold on
  shadedErrorBar(1:25,nanmean(r_df,2),nanstd(r_df,[],2)/sqrt(50),'r');
  axis([1 25 -0.35 0.35]); tp_editplots
  line([1 25],[0 0],'color',[.5 .5 .5],'linestyle',':')
  h=p<(fdr1(p,0.05,0)/54);
  plot(find(h),mean(r(find(h),:),2),'k.','markersize',8)
  h=p<(fdr1(p_df,0.05,0)/54);
  if ~isnan(nanmean(r_df(find(h),:),2))
  plot(find(h),nanmean(r_df(find(h),:),2),'k.','markersize',8)
  end
  axis square
  set(gca,'tickdir','out','xtick',1:4:25,'xticklabel',freqoi(1:4:25))
  set(gca,'tickdir','out','ytick',[-0.35 0 0.35],'yticklabel',[-0.35 0 0.35])

end
set(gcf,'Position',[50 50  800 1200])
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_receptors_abs%d_dd%s_nr2_v%d.pdf',isabs,num2str(dd),v))
close

%% RECEPTORS: DERIVATIVE OF PUPIL 

orig = [46 64 37];
locs_mni=2*(locs-repmat(orig,[size(locs,1) 1]))/10;
isabs=1;
  
dd = 1;
for irecep = 1 : 54

  freqoi=2.^(1:(1/4):7);
  par = nan(8799,1);
  par(ismember(BNA.tissue_5mm,[1:2:246]))=spatfiltergauss(nanmean(gdat(:,irecep),2),locs_mni,dd,BNA.grid_5mm(ismember(BNA.tissue_5mm,[1:2:246]),:)/10);
  clear r p par2 par1
  for isubj = 1 : 50
    isubj
    k=0;
    pp = plt_all.corr_src_df(:,:,isubj);
    for igrid = 1:2:246
      k=k+1;
      par1(k) = mean(par(BNA.tissue_5mm == igrid));
    end
    k = 0;
    for igrid = 1 : 246
      k= k+1;
      tmp(k,:) = mean(pp(BNA.tissue_5mm == igrid,:));
    end
    par2=(tmp(1:2:246,:)+tmp(2:2:246,:))/2;
    %     par2(par2<0)=0;
    if isabs==1
      [r(:,isubj),p(:,isubj)] = corr(par1(:),abs(par2),'type','pearson');
    else
      [r(:,isubj),p(:,isubj)] = corr(par1(:),par2,'type','pearson');
    end
  end
  

  [~,p]=ttest(r,zeros(size(r)),'dim',2);

  if irecep == 1
    figure_w; 
    subplot(6,5,irecep); hold on
  elseif irecep ~= 1 && irecep<28
    subplot(6,5,irecep); hold on
  end
  
  if irecep == 28
    set(gcf,'Position',[50 50 800 1200])
    print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_receptors_abs%d_dd%s_nr1_v%d.pdf',isabs,num2str(dd),v))
    close
    figure_w; 
    subplot(6,5,irecep-27); hold on
  elseif irecep ~= 28 && irecep>27
    subplot(6,5,irecep-27); hold on
  end
    
  
  title(sprintf('%s',lab{irecep}))

  shadedErrorBar(1:25,mean(r,2),std(r,[],2)/sqrt(50),'k');
  axis([1 25 -0.35 0.35]); tp_editplots
  line([1 25],[0 0],'color',[.5 .5 .5],'linestyle',':')
  h=p<(fdr1(p,0.05,0)/54);
  plot(find(h),mean(r(find(h),:),2),'k.','markersize',8)
  axis square
  set(gca,'tickdir','out','xtick',1:4:25,'xticklabel',freqoi(1:4:25))
  set(gca,'tickdir','out','ytick',[-0.35 0 0.35],'yticklabel',[-0.35 0 0.35])

end
set(gcf,'Position',[50 50  800 1200])
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_receptors_derivative_abs%d_dd%s_nr2_v%d.pdf',isabs,num2str(dd),v))
close

%% PLOT CROSS CORRELATION
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cols = cbrewer('seq', 'Oranges', 45,'pchip');
cols = cols(1:end-15,:); cols=cols(end:-1:1,:);

v=1
SUBJLIST=1:41; SUBJLIST([10,12,17,19,22,27,35,38,39,40])=[];
i = 0;
for isubj = SUBJLIST
    i = i + 1; i
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_mue_src_pupil_power_correlations_s%d_b1_v%d.mat',isubj,v));
    
    for ifreq = 1 : 25
        mue.xcorr{ifreq}(:,i) =  nanmean(outp.xcorr{ifreq},2);
        mue.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr{ifreq},2);
        
    end
end

% plt_mue.xcorr_ord= zeros(size(ser_hh,2)-1,25,size(plt_mue.corr_sens,3));
% for i = 1 : size(ser_hh,2)-1
%   idx = lay.pos(:,2)>ser_hh(i) & lay.pos(:,2)<ser_hh(i+1);
%   for ifreq = 1 : 25
%   mue.xcorr_ord{ifreq}(i,:) = (mue.xcorr{ifreq}(idx,:),1);
%   end
% end

figure; set(gcf,'color','w') ;
subplot(2,2,3); hold on; title('Muenster')

for ifreq = 1 : 25
    plot(outp.xcorr_lags{ifreq},nanmean(mue.xcorr{ifreq},2),'color',cols(ifreq,:))
end

axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

cols = cbrewer('seq', 'Reds', 35,'pchip');
cols = cols(1:end-5,:); cols=cols(end:-1:1,:);
v=1; SUBJLIST = 1:24; SUBJLIST([5,9])=[];
i = 0;
for isubj = SUBJLIST
    i = i + 1;
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_gla_src_pupil_power_correlations_s%d_b1_v%d.mat',isubj,v)); 
    for ifreq = 1 : 25
        gla.xcorr{ifreq}(:,i) =  nanmean(outp.xcorr{ifreq},2);
        gla.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr_df{ifreq},2);
    end
end

subplot(2,2,1); hold on; title('Glasgow')

for ifreq = 1 : 25
    plot(outp.xcorr_lags{ifreq},nanmean(gla.xcorr{ifreq},2),'color',cols(ifreq,:))
end

axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

cols = cbrewer('seq', 'Blues', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
v=1
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
i = 0;
for isubj = SUBJLIST
    i = i + 1;
    d = dir(sprintf('/home/tpfeffer/pp/proc/src/pp_src_pupil_power_correlations_s%d_b*_v%d.mat',isubj,v));
    if length(d)==1
        load([d(1).folder '/' d(1).name])   
        for ifreq = 1 : 25       
            hh.xcorr{ifreq}(:,i)    =  nanmean(outp.xcorr{ifreq},2);
            hh.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr_df{ifreq},2);
        end
    else
        load([d(1).folder '/' d(1).name])
        for ifreq = 1 : 25
            hh.xcorr{ifreq}(:,i)    =  nanmean(outp.xcorr{ifreq},2)./2;
            hh.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr_df{ifreq},2)./2;
        end
        load([d(2).folder '/' d(2).name])
        for ifreq = 1 : 25
            hh.xcorr{ifreq}(:,i)    = hh.xcorr{ifreq}(:,i)+ nanmean(outp.xcorr{ifreq},2)./2;
            hh.xcorr_df{ifreq}(:,i) = hh.xcorr_df{ifreq}(:,i)+nanmean(outp.xcorr_df{ifreq},2)./2;
        end
    end
end
% figure; set(gcf,'color','w') ;

subplot(2,2,2); hold on; title('Hamburg')

for ifreq = 1 : 25
    plot(outp.xcorr_lags{ifreq},nanmean(hh.xcorr{ifreq},2),'color',cols(ifreq,:))
end
axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

% cfg=[];
% cfg.layout='CTF275.lay';
% lay = ft_prepare_layout(cfg);
% % lay.pos([203 276 277],:)=[];
% minmax_hh=[min(lay.pos(:,2)) max(lay.pos(:,2))];
% ser_hh = linspace(minmax_hh(1),minmax_hh(2),40);
% 
% % for ifreq = 1 : 25
% %     hh.xcorr_ord{ifreq} = nan(size(hh.xcorr{ifreq},1),size(hh.xcorr{ifreq},2));
% for i = 3 : size(ser_hh,2)-1
%   idx = lay.pos(:,2)>ser_hh(i) & lay.pos(:,2)<ser_hh(i+1);
%   for ifreq = 1 : 25
%        xcorr_ord{ifreq}(i,:,:) = squeeze(nanmean(hh.xcorr{ifreq}(:,idx,:),2));
%         [~,j]=max(nanmean(nanmean(squeeze(xcorr_ord{ifreq}(i,:,:)),2),3));
%         lag_max(i,ifreq)=outp.xcorr_lags{ifreq}(j);
%         [~,j]=min(nanmean(nanmean(squeeze(xcorr_ord{ifreq}(i,:,:)),2),3));
%         lag_min(i,ifreq)=outp.xcorr_lags{ifreq}(j);
%   end
% end

% cfg=[];
% cfg.layout='4D248.lay';
% lay = ft_prepare_layout(cfg);
% minmax_gla=[min(lay.pos(1:248,2)) max(lay.pos(1:248,2))];
% ser_gla = linspace(minmax_gla(1),minmax_gla(2),40);
% plt_gla.corr_sens_ord= zeros(size(ser_gla,2)-1,25,22);
% 
% for i = 1 : size(ser_gla,2)-1
%   idx = lay.pos(:,2)<ser_gla(i+1) & lay.pos(:,2)>ser_gla(i);
%   plt_gla.corr_sens_ord(i,:,:) = nanmean(plt_gla.corr_sens(idx,:,:),1);
% end
% 
% cfg=[];
% cfg.layout='CTF275.lay';
% lay = ft_prepare_layout(cfg);
% % lay.pos([203 276 277],:)=[];
% minmax_hh=[min(lay.pos(:,2)) max(lay.pos(:,2))];
% ser_hh = linspace(minmax_hh(1),minmax_hh(2),40);
% 
% plt_mue.corr_sens_ord= zeros(size(ser_hh,2)-1,25,size(plt_mue.corr_sens,3));
% for i = 1 : size(ser_hh,2)-1
%   idx = lay.pos(:,2)>ser_hh(i) & lay.pos(:,2)<ser_hh(i+1);
%   plt_mue.corr_sens_ord(i,:,:) = nanmean(plt_mue.corr_sens(idx,:,:),1);
% end
% plt_hh.corr_sens_ord= zeros(size(ser_hh,2)-1,25,28);
% for i = 1 : size(ser_hh,2)-1
%   idx = lay.pos(:,2)>ser_hh(i) & lay.pos(:,2)<ser_hh(i+1);
%   plt_hh.corr_sens_ord(i,:,:) = nanmean(plt_hh.corr_sens(idx,:,:),1);
% end
% 


cols = cbrewer('seq', 'Greys', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
for ifreq = 1 : 25
    
    all.xcorr{ifreq}=cat(2,gla.xcorr{ifreq},mue.xcorr{ifreq},hh.xcorr{ifreq});
     all.xcorr_df{ifreq}=cat(2,gla.xcorr_df{ifreq},mue.xcorr_df{ifreq},hh.xcorr_df{ifreq});
   
end


subplot(2,2,4); hold on; title('Pooled')

for ifreq = 1 : 25
plot(outp.xcorr_lags{ifreq},nanmean(all.xcorr{ifreq},2),'color',cols(ifreq,:))

end

axis([-5 5 -0.06 0.04])
xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

%% XCORR FOR SELECTED FREQUENCIES
i = 0;  clear sig
for ifreq = 21:24
    i = i + 1;
    sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{25},1)));
end

sig = zscore(nanmean(nanmean(sig,2),3));
[i,j]= max(sig); lag = outp.xcorr_lags{25}(j);

figure_w
subplot(2,2,1); hold on
plot(outp.xcorr_lags{25},sig,'color',cols(25,:))
line([lag lag],[0 i],'color',cols(25,:))
lag
[i,j]= max(sig(1:find(outp.xcorr_lags{25}<0.85,1,'last'))); lag = outp.xcorr_lags{25}(j);
line([lag lag],[0 i],'color',cols(25,:))
lag
clear sig
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];

i = 0;
for ifreq = 1:4
    i = i + 1;
    sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{5},1)));
end
    
sig = zscore(nanmean(nanmean(sig,2),3));
[i,j]= min(sig); lag = outp.xcorr_lags{5}(j);
lag
plot(outp.xcorr_lags{5},sig,'color',cols(1,:))
line([lag lag],[i 0],'color',cols(1,:))

clear sig

i = 0; 
for ifreq = 11:13
    i = i + 1;
    sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{14},1)));
end
sig = zscore(nanmean(nanmean(sig,2),3));
[i,j]= max(sig); lag = outp.xcorr_lags{14}(j);

plot(outp.xcorr_lags{14},sig,'color',cols(14,:))
axis([-5 5 -4 4]);xlabel('Lag [s]'); ylabel('Correlation coeff. (z-scored)');
tp_editplots;
line([lag lag],[0 i],'color',cols(14,:))
lag
line([0 0],[-4.5 4.5],'color',[.8 .8 .8],'linestyle',':')
line([-5 5],[0 0],'color',[0.8 .8 .8],'linestyle',':')
% cfg=[];
    
    
% SAME FOR DERIVATIVE
i = 0;  clear sig
for ifreq = 21:24
    i = i + 1;
    sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{25},1)));
end

sig = zscore(nanmean(nanmean(sig,2),3));
    
subplot(2,2,2); hold on
plot(outp.xcorr_lags{25},sig,'color',cols(25,:))
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];

clear sig

i = 0;
for ifreq = 1:4
    i = i + 1;
    sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{5},1)));
end
    
sig = zscore(nanmean(nanmean(sig,2),3));

plot(outp.xcorr_lags{5},sig,'color',cols(1,:))

clear sig

i = 0; 
for ifreq = 11:13
    i = i + 1;
    sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{14},1)));
end
sig = zscore(nanmean(nanmean(sig,2),3));
    
plot(outp.xcorr_lags{14},sig,'color',cols(14,:))
axis([-5 5 -4 4]);xlabel('Lag [s]'); ylabel('Correlation coeff. (z-scored)');
tp_editplots; 
line([0 0],[-4 4],'color',[.8 .8 .8],'linestyle',':')
line([-5 5],[0 0],'color',[0.8 .8 .8],'linestyle',':')
% cfg=[];
colorbar
    
    



