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
v = 2

[plt_gla,plt_hh,plt_mue,plt_hh_cnt,plt_all]=pp_load_results(v);

colors = cbrewer('qual', 'Set3', 10,'pchip'); 
colors = colors(4:6,:);

%% PLOT MUTUAL INFORMATION
freqoi=2.^(1:(1/4):7); 

figure_w

subplot(2,3,1);
imagesc(nanmean(plt_gla.sens_mi_ord,3),[0 .005])
colormap(plasma); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,2);
imagesc(nanmean(plt_hh.sens_mi_ord,3),[0 .005])
colormap(plasma); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,3);
imagesc(nanmean(plt_mue.sens_mi_ord,3),[0 .005])
colormap(plasma); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,4);
par = cat(3,plt_gla.sens_mi_ord,plt_hh.sens_mi_ord,plt_mue.sens_mi_ord);
imagesc(nanmean(plt_mue.sens_mi_ord,3),[0 .005])
colormap(plasma); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,5);
imagesc(nanmean(plt_hh_cnt.sens_mi_ord,3),[0 .005])
colormap(plasma); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_MI_v%d.pdf',v))

freqoi=2.^(1:(1/4):7); 

tot_size = size(plt_mue.sens_mi,3)+size(plt_gla.sens_mi,3)+size(plt_hh.sens_mi,3);

par_gla = mean(nanmean(plt_gla.sens_mi,1),3);
par_hh  = mean(nanmean(plt_hh.sens_mi,1),3);
par_mue = nanmean(nanmean(plt_mue.sens_mi,1),3);
par_hh_cnt  = mean(nanmean(plt_hh_cnt.sens_mi,1),3);

std_gla = std(nanmean(plt_gla.sens_mi,1),[],3)/sqrt(size(plt_gla.sens_mi,3));
std_hh  = std(nanmean(plt_hh.sens_mi,1),[],3)/sqrt(size(plt_hh.sens_mi,3));
std_mue = nanstd(nanmean(plt_mue.sens_mi,1),[],3)/sqrt(size(plt_mue.sens_mi,3));
std_hh_cnt  = std(nanmean(plt_hh_cnt.sens_mi,1),[],3)/sqrt(size(plt_hh_cnt.sens_mi,3));

par_all = mean(cat(2,squeeze(nanmean(plt_gla.sens_mi,1)),squeeze(nanmean(plt_hh.sens_mi,1)),squeeze(nanmean(plt_mue.sens_mi,1))),2);
std_all = std(cat(2,squeeze(nanmean(plt_gla.sens_mi,1)),squeeze(nanmean(plt_hh.sens_mi,1)),squeeze(nanmean(plt_mue.sens_mi,1))),[],2)/sqrt(tot_size);

figure_w
subplot(4,3,1)
shadedErrorBar(log10(freqoi),par_gla,std_gla,{'color',colors(1,:)})
axis([.3 2.11 0 0.007])
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,2); hold on 
shadedErrorBar(log10(freqoi),par_hh,std_hh,{'color',colors(2,:)})
axis([.3 2.11 0 0.0030])
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

shadedErrorBar(log10(freqoi),par_hh_cnt,std_hh_cnt,{'color',colors(2,:),'linestyle',':'})
axis([.3 2.11 0 0.0030])
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

% stats task vs. rest
[~,p]=ttest(nanmean(plt_hh.sens_mi,1),nanmean(plt_hh_cnt.sens_mi,1),'dim',3); h = p< fdr1(p(:),0.05,0);
plot(log10(freqoi(h)),0.0025*ones(sum(h),1),'r*','markersize',2)
plot(log10(freqoi(p<0.05)),0.003*ones(sum(p<0.05),1),'k*','markersize',2)

subplot(4,3,3)
shadedErrorBar(log10(freqoi),par_mue,std_mue,{'color',colors(3,:)})
axis([.3 2.11 0 0.0060])
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,4)
shadedErrorBar(log10(freqoi),par_all,std_all,{'color','k'})
axis([.3 2.11 0 0.0050])
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))


print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_MI_lineplot_v%d.pdf',v))
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

subplot(2,3,4);
imagesc(nanmean(plt_hh_cnt.corr_sens_ord(3:end,:,:),3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

pooled = cat(3,plt_hh.corr_sens_ord,plt_gla.corr_sens_ord,plt_mue.corr_sens_ord);

subplot(2,3,5);
imagesc(nanmean(pooled(3:end,:,:),3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3);
h = p<fdr1(p(:),0.05,0);

subplot(2,3,6);
imagesc(nanmean(pooled(3:end,:,:),3).*h(3:end,:),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_anterior_post_v%d.pdf',v))

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

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

figure_w
subplot(2,4,1)
[h,p]=ttest(plt_gla.corr_src_BNA(idx_sorted,:,:),zeros(size(plt_gla.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_gla.corr_src_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,2)
[h,p]=ttest(plt_hh.corr_src_BNA(idx_sorted,:,:),zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_hh.corr_src_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,3)
[h,p]=ttest(plt_mue.corr_src_BNA(idx_sorted,:,:),zeros(size(plt_mue.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_mue.corr_src_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,4)

pooled= cat(3,plt_hh.corr_src_BNA(idx_sorted,:,:),plt_gla.corr_src_BNA(idx_sorted,:,:),plt_mue.corr_src_BNA(idx_sorted,:,:)); 

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,0);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,5)
[h,p]=ttest(plt_gla.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_gla.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_gla.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,6)
[h,p]=ttest(plt_hh.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_hh.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,7)
[h,p]=ttest(plt_mue.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_mue.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_mue.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,8)

pooled= cat(3,plt_hh.corr_src_df_BNA(idx_sorted,:,:),plt_gla.corr_src_df_BNA(idx_sorted,:,:),plt_mue.corr_src_df_BNA(idx_sorted,:,:)); 

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
imagesc(nanmean(plt_hh.corr_src_BNA_cnt(idx_sorted,:,:),3),[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,2)
imagesc(nanmean(plt_hh.corr_src_BNA_cnt(idx_sorted,:,:),3)-nanmean(plt_hh.corr_src_BNA(idx_sorted,:,:),3),[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,3)
[h,p]=ttest(plt_hh.corr_src_BNA_cnt,zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc((nanmean(plt_hh.corr_src_BNA_cnt(idx_sorted,:,:),3)-nanmean(plt_hh.corr_src_BNA(idx_sorted,:,:),3)).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,4)
[h,p]=ttest(plt_hh.corr_src_BNA_cnt(idx_sorted,:,:),plt_hh.corr_src_BNA(idx_sorted,:,:),'dim',3); 
h=p<0.05;
imagesc((nanmean(plt_hh.corr_src_BNA_cnt(idx_sorted,:,:),3)-nanmean(plt_hh.corr_src_BNA(idx_sorted,:,:),3)).*h,[-0.05 0.05])
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
    mue.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr_df{ifreq},2);
  end
end

figure; set(gcf,'color','w') ;
subplot(2,2,3); hold on; title('Muenster')

for ifreq = 1 : 25
  plot(outp.xcorr_lags{ifreq},nanmean(mue.xcorr_df{ifreq},2),'color',cols(ifreq,:))
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
  plot(outp.xcorr_lags{ifreq},nanmean(gla.xcorr_df{ifreq},2),'color',cols(ifreq,:))
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
  plot(outp.xcorr_lags{ifreq},nanmean(hh.xcorr_df{ifreq},2),'color',cols(ifreq,:))
end
axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

cols = cbrewer('seq', 'Greys', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
for ifreq = 1 : 25
  all.xcorr{ifreq}=cat(2,gla.xcorr{ifreq},mue.xcorr{ifreq},hh.xcorr{ifreq});
  all.xcorr_df{ifreq}=cat(2,gla.xcorr_df{ifreq},mue.xcorr_df{ifreq},hh.xcorr_df{ifreq});
end

subplot(2,2,4); hold on; title('Pooled')

for ifreq = 1 : 25
  plot(outp.xcorr_lags{ifreq},nanmean(all.xcorr_df{ifreq},2),'color',cols(ifreq,:))
  
end

axis([-5 5 -0.06 0.04])
xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_df_allfreqs_v%d.pdf',v))

%% XCORR FOR SELECTED FREQUENCIES
task = 0;

if task 
cols = cbrewer('seq', 'Blues', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);    
else
cols = cbrewer('seq', 'Greys', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);    
end

figure_w
subplot(2,2,1); hold on

i = 0; clear sig
for ifreq = 1:4
    i = i + 1;
       if task
        sig(:,:,i) = interp1(1:size(hh_cnt.xcorr{ifreq},1),hh_cnt.xcorr{ifreq},linspace(1,size(hh_cnt.xcorr{ifreq},1),size(hh_cnt.xcorr{5},1)));
       else
        sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{5},1)));
   end
end
    
[c,p]=permutest(nanmean(sig,3),zeros(size(nanmean(sig,3))),1,0.05,1000,1);

sig = zscore(nanmean(nanmean(sig,2),3));
[i,j]= min(sig); lag = outp.xcorr_lags{5}(j);
lag
plot(outp.xcorr_lags{5},sig,'color',cols(1,:))

for i = 1 : length(c)
    if p(i)>0.05; continue; end
    plot(outp.xcorr_lags{5}(c{i}),sig(c{i}),'*','markersize',3,'color',cols(5,:));
end
line([lag lag],[i 0],'color',cols(1,:))

i = 0; clear sig
for ifreq = 11:13
    i = i + 1;
    if task
        sig(:,:,i) = interp1(1:size(hh_cnt.xcorr{ifreq},1),hh_cnt.xcorr{ifreq},linspace(1,size(hh_cnt.xcorr{ifreq},1),size(hh_cnt.xcorr{14},1)));
    else
        sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{14},1)));
    end
end

[c,p]=permutest(nanmean(sig,3),zeros(size(nanmean(sig,3))),1,0.05,1000,1);

sig = zscore(nanmean(nanmean(sig,2),3));
[i,j]= max(sig); lag = outp.xcorr_lags{14}(j);

plot(outp.xcorr_lags{14},sig,'color',cols(14,:))
for i = 1 : length(c)
    if p(i)>0.05; continue; end
    plot(outp.xcorr_lags{14}(c{i}),sig(c{i}),'*','markersize',3,'color',cols(14,:));
end
axis([-5 5 -6 4]);xlabel('Lag [s]'); ylabel('Correlation coeff. (z-scored)');
tp_editplots;
line([lag lag],[0 i],'color',cols(14,:))
lag
line([0 0],[-6 4],'color',[.8 .8 .8],'linestyle',':')
line([-5 5],[0 0],'color',[0.8 .8 .8],'linestyle',':')

i = 0; clear sig
for ifreq = 21:24
    i = i + 1;
    if task 
        sig(:,:,i) = interp1(1:size(hh_cnt.xcorr{ifreq},1),hh_cnt.xcorr{ifreq},linspace(1,size(hh_cnt.xcorr{ifreq},1),size(hh_cnt.xcorr{25},1)));
    else
        sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{25},1)));
    end
end

[c,p]=permutest(nanmean(sig,3),zeros(size(nanmean(sig,3))),1,0.05,1000,1);

sig = zscore(nanmean(nanmean(sig,2),3));
[i,j]= max(sig); lag = outp.xcorr_lags{25}(j);
plot(outp.xcorr_lags{25},sig,'color',cols(25,:))
for i = 1 : length(c)
    if p(i)>0.05; continue; end
    plot(outp.xcorr_lags{25}(c{i}),sig(c{i}),'*','markersize',3,'color',cols(25,:));
end
line([lag lag],[0 i],'color',cols(25,:))
lag
[i,j]= max(sig(1:find(outp.xcorr_lags{25}<0.85,1,'last'))); lag = outp.xcorr_lags{25}(j);
line([lag lag],[0 i],'color',cols(25,:))
lag
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];

% SAME FOR DERIVATIVE
i = 0;  clear sig
for ifreq = 21:24
    i = i + 1;
    if task
        sig(:,:,i) = interp1(1:size(hh_cnt.xcorr_df{ifreq},1),hh_cnt.xcorr_df{ifreq},linspace(1,size(hh_cnt.xcorr_df{ifreq},1),size(hh_cnt.xcorr_df{25},1)));
    else
        sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{25},1)));
   
    end
end

[c,p]=permutest(nanmean(sig,3),zeros(size(nanmean(sig,3))),1,0.05,1000,1);

sig = zscore(nanmean(nanmean(sig,2),3));
    
subplot(2,2,2); hold on
plot(outp.xcorr_lags{25},sig,'color',cols(25,:))
for i = 1 : length(c)
    if p(i)>0.05; continue; end
    plot(outp.xcorr_lags{25}(c{i}),sig(c{i}),'*','markersize',3,'color',cols(25,:));
end
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];

clear sig

i = 0;
for ifreq = 1:4
    i = i + 1;
    if task
        sig(:,:,i) = interp1(1:size(hh_cnt.xcorr_df{ifreq},1),hh_cnt.xcorr_df{ifreq},linspace(1,size(hh_cnt.xcorr_df{ifreq},1),size(hh_cnt.xcorr_df{5},1)));
    else
    sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{5},1)));
    end
end

[c,p]=permutest(nanmean(sig,3),zeros(size(nanmean(sig,3))),1,0.05,1000,1);

sig = zscore(nanmean(nanmean(sig,2),3));

plot(outp.xcorr_lags{5},sig,'color',cols(1,:))
for i = 1 : length(c)
    if p(i)>0.05; continue; end
    plot(outp.xcorr_lags{5}(c{i}),sig(c{i}),'*','markersize',3,'color',cols(5,:));
end
clear sig

i = 0; 
for ifreq = 11:13
    i = i + 1;
    if task
        sig(:,:,i) = interp1(1:size(hh_cnt.xcorr_df{ifreq},1),hh_cnt.xcorr_df{ifreq},linspace(1,size(hh_cnt.xcorr_df{ifreq},1),size(hh_cnt.xcorr_df{14},1)));
    else
    sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{14},1)));
    end
end
[c,p]=permutest(nanmean(sig,3),zeros(size(nanmean(sig,3))),1,0.05,1000,1);

sig = zscore(nanmean(nanmean(sig,2),3));
    
plot(outp.xcorr_lags{14},sig,'color',cols(14,:))
for i = 1 : length(c)
    if p(i)>0.05; continue; end
    plot(outp.xcorr_lags{14}(c{i}),sig(c{i}),'*','markersize',3,'color',cols(14,:));
end
axis([-5 5 -6 4]);xlabel('Lag [s]'); ylabel('Correlation coeff. (z-scored)');
tp_editplots; 
line([0 0],[-6 4],'color',[.8 .8 .8],'linestyle',':')
line([-5 5],[0 0],'color',[0.8 .8 .8],'linestyle',':')
% cfg=[];
colorbar
    
if task
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_cnt_xcorr_pooledfreqs_v%d.pdf',v))
else
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_pooledfreqs_v%d.pdf',v))
end
%%

cols = cbrewer('seq', 'Blues', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
v=1
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
i = 0;
for isubj = SUBJLIST
    i = i + 1;
    d = dir(sprintf('/home/tpfeffer/pp/proc/src/pp_cnt_src_pupil_power_correlations_s%d_b*_v%d.mat',isubj,v));
    if length(d)==1
        load([d(1).folder '/' d(1).name])   
        for ifreq = 1 : 25       
            hh_cnt.xcorr{ifreq}(:,i)    =  nanmean(outp.xcorr{ifreq},2);
            hh_cnt.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr_df{ifreq},2);
        end
    else
        load([d(1).folder '/' d(1).name])
        for ifreq = 1 : 25
            hh_cnt.xcorr{ifreq}(:,i)    =  nanmean(outp.xcorr{ifreq},2)./2;
            hh_cnt.xcorr_df{ifreq}(:,i) = nanmean(outp.xcorr_df{ifreq},2)./2;
        end
        load([d(2).folder '/' d(2).name])
        for ifreq = 1 : 25
            hh_cnt.xcorr{ifreq}(:,i)    = hh_cnt.xcorr{ifreq}(:,i)+ nanmean(outp.xcorr{ifreq},2)./2;
            hh_cnt.xcorr_df{ifreq}(:,i) = hh_cnt.xcorr_df{ifreq}(:,i)+nanmean(outp.xcorr_df{ifreq},2)./2;
        end
    end
end

figure; set(gcf,'color','w') ;
subplot(2,2,1); hold on; title('Hamburg')
for ifreq = 1 : 25
    plot(outp.xcorr_lags{ifreq},nanmean(hh_cnt.xcorr_df{ifreq},2),'color',cols(ifreq,:))
end
axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([0.93 0.93],[-0.06 0.04],'color','k')
line([0.52 0.52],[-0.06 0.04],'color','k')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_cnt_xcorr_df_allfreqs_v%d.pdf',v))

%% PLOT SCALING EXPONENTS

% load fooof results and power spectra
% fooof = pp_load_fooof_results(v);
% fooof.ps_* = power spectra
% fooof.psfit_* = fooof fits of spectra
% fooof.aper_* = offset and slope of 1/f component

% ----------------------------------------
% PLOT OFFSET AND SLOPE OF FOOOF FITS
% ----------------------------------------

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

figure_w;

subplot(4,3,1); hold on; title('Offset - HH')
par = squeeze(nanmean(fooof.aper_hh(1,idx_sorted,3,:)-fooof.aper_hh(1,idx_sorted,1,:),3));
plot(nanmean(par,2),'k'); 
d = squeeze(fooof.aper_hh(1,idx_sorted,3,:)-fooof.aper_hh(1,idx_sorted,1,:));
[h,p] = ttest(d,zeros(size(d)),'dim',2); 
plot(find(h),nanmean(par(find(h),:),2),'r.','markersize',5)
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in offs.'); xlabel('Region')

subplot(4,3,4); hold on; title('Slope - HH')
par = squeeze(nanmean(fooof.aper_hh(2,idx_sorted,3,:)-fooof.aper_hh(2,idx_sorted,1,:),3));
plot(nanmean(par,2),'k'); 
d = squeeze(fooof.aper_hh(2,idx_sorted,3,:)-fooof.aper_hh(2,idx_sorted,1,:));
[h,p] = ttest(d,zeros(size(d)),'dim',2); 
plot(find(h),nanmean(par(find(h),:),2),'r.','markersize',5)
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in slope'); xlabel('Region')

subplot(4,3,2); hold on; title('Offset - GLA')
par = squeeze(nanmean(fooof.aper_gla(1,idx_sorted,3,:)-fooof.aper_gla(1,idx_sorted,1,:),3));
plot(nanmean(par,2),'k'); 
d = squeeze(fooof.aper_gla(1,idx_sorted,3,:)-fooof.aper_gla(1,idx_sorted,1,:));
[h,p] = ttest(d,zeros(size(d)),'dim',2); 
plot(find(h),nanmean(par(find(h),:),2),'r.','markersize',5)
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in offs.'); xlabel('Region')

subplot(4,3,5); hold on; title('Slope - GLA')
par = squeeze(nanmean(fooof.aper_gla(2,idx_sorted,3,:)-fooof.aper_gla(2,idx_sorted,1,:),3));
plot(nanmean(par,2),'k'); 
d = squeeze(fooof.aper_gla(2,:,3,:)-fooof.aper_gla(2,:,1,:));
[h,p] = ttest(d,zeros(size(d)),'dim',2); 
plot(find(h),nanmean(par(find(h),:),2),'r.','markersize',5)
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in slope'); xlabel('Region')

subplot(4,3,3); hold on; title('Offset - MUE')
par = squeeze(nanmean(fooof.aper_mue(1,idx_sorted,3,:)-fooof.aper_mue(1,idx_sorted,1,:),3));
plot(nanmean(par,2),'k'); 
d = squeeze(fooof.aper_mue(1,idx_sorted,3,:)-fooof.aper_mue(1,idx_sorted,1,:));
[h,p] = ttest(d,zeros(size(d)),'dim',2); 
plot(find(h),nanmean(par(find(h),:),2),'r.','markersize',5)
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in offs.'); xlabel('Region')

subplot(4,3,6); hold on; title('Slope - MUE')
par = squeeze(nanmean(fooof.aper_mue(2,idx_sorted,3,:)-fooof.aper_mue(2,idx_sorted,1,:),3));
plot(nanmean(par,2),'k'); 
d = squeeze(fooof.aper_mue(2,idx_sorted,3,:)-fooof.aper_mue(2,idx_sorted,1,:));
[h,p] = ttest(d,zeros(size(d)),'dim',2); 
plot(find(h),nanmean(par(find(h),:),2),'r.','markersize',5)
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in slope'); xlabel('Region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_fooof_offset_and_slope_v%d.pdf',v))

%% CORRELATIONS
d = pow_hi_hh-pow_lo_hh;
idx = ~isnan(pow_lo_hh(:,1,1));
for isubj = 1 : 28
  tmp1=d(idx,:,isubj);
  tmp2=d(idx,:,isubj);

  r(:,:,isubj)=corr(tmp1',tmp2');
  
end

h=ttest(r,zeros(size(r)),'dim',3);

figure_w; 
subplot(3,2,1);
imagesc(nanmean(r,3),[-0.5 0.5]); colormap(cmap); axis square
subplot(3,2,2);
imagesc(nanmean(r,3).*h,[-0.5 0.5]); colormap(cmap); axis  square

d = pow_hi_gla-pow_lo_gla;
idx = ~isnan(pow_lo_hh(:,1,1));
for isubj = 1 : 22
  tmp1=d(idx,:,isubj);
  tmp2=d(idx,:,isubj);

  r(:,:,isubj)=corr(tmp1',tmp2');
  
end

h=ttest(r,zeros(size(r)),'dim',3);

subplot(3,2,3);
imagesc(nanmean(r,3),[-0.5 0.5]); colormap(cmap); axis square
subplot(3,2,4);
imagesc(nanmean(r,3).*h,[-0.5 0.5]); colormap(cmap); axis  square

d1 = pow_hi_hh-pow_lo_hh;
d2 = pow_hi_gla-pow_lo_gla;
d = cat(3,d1,d2);

idx = ~isnan(pow_lo_hh(:,1,1));
for isubj = 1 : 50
  tmp1=d(idx,:,isubj);
  tmp2=d(idx,:,isubj);

  r(:,:,isubj)=corr(tmp1',tmp2');
  
end

[h,~,~,s]=ttest(r,zeros(size(r)),'dim',3);

subplot(3,2,5);
imagesc(nanmean(r,3),[-0.5 0.5]); colormap(cmap); axis square
subplot(3,2,6);
imagesc(nanmean(r,3).*h,[-0.5 0.5]); colormap(cmap); axis  square

ffx = 2:0.5:128;
ffx(ffx>48 & ffx<52) = [];
ffx(ffx>98 & ffx<102) = [];
figure_w; hold on

plot(log10(ffx),100*sum(h>0&s.tstat<0)/239)
plot(log10(ffx),100*sum(h>0&s.tstat>0)/239)
tp_editplots

set(gca,'xtick',log10(ffx([1,5,13,29,61,118,239])),'xticklabel',ffx([1,5,13,29,61,118,239]))
xlabel('Frequency [Hz]'); ylabel ('Fraction of sign. correlations')
