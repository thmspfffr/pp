%%

clear
restoredefaultpath

addpath ~/pconn/matlab/
addpath ~/pp/matlab/
load /home/gnolte/meth/templates/mri.mat
load /home/gnolte/meth/templates/sa_template.mat
addpath ~/Documents/MATLAB/fieldtrip-20190224/
ft_defaults



addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
v = 2;

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
imagesc(nanmean(plt_gla.corr_sens_ord,3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,2);
imagesc(nanmean(plt_hh.corr_sens_ord,3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,3);
imagesc(nanmean(plt_mue.corr_sens_ord,3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,4);
imagesc(nanmean(plt_hh_cnt.corr_sens_ord,3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

pooled = cat(3,plt_hh.corr_sens_ord,plt_gla.corr_sens_ord,plt_mue.corr_sens_ord);

subplot(2,3,5);
imagesc(nanmean(pooled,3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3);
h = p<fdr1(p(:),0.05,1);
subplot(2,3,6);
imagesc(nanmean(pooled,3).*h,[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_anterior_post_v%d.pdf',v))

figure_w;
% TASK VS REST

[h,p]=ttest(plt_hh_cnt.corr_sens_ord,zeros(size(plt_hh_cnt.corr_sens_ord)),'dim',3);
h = p<fdr1(p(:),0.05,1); 
% h = p < 0.01;
par = nanmean(plt_hh_cnt.corr_sens_ord,3);
subplot(2,3,4);
imagesc(par.*h,[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

[h,p]=ttest(plt_hh_cnt.corr_sens_ord,plt_hh.corr_sens_ord,'dim',3);
h = p<fdr1(p(:),0.05,1); 
h = p < 0.01;
par = nanmean(plt_hh_cnt.corr_sens_ord-plt_hh.corr_sens_ord,3);
subplot(2,3,5);
imagesc(par,[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

subplot(2,3,6);
imagesc(par.*h,[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_task_sens_anterior_post_v%d.pdf',v))


%% PLOT AVERAGES ACROSS SPACE, HAMBURG, GLASGOW, ALL

% -----------------
% SENSOR SPACE - POWER
% -----------------
freqoi=2.^(1:(1/4):7);

tot_size = size(plt_mue.corr_sens,3)+size(plt_gla.corr_sens,3)+size(plt_hh.corr_sens,3);

par_gla = mean(log10(nanmean(plt_gla.pow_sens,1)),3);
par_hh  = mean(log10(nanmean(plt_hh.pow_sens,1)),3);
par_mue = nanmean(log10(nanmean(plt_mue.pow_sens,1)),3);
par_cnt = nanmean(log10(nanmean(plt_hh_cnt.pow_sens,1)),3);

std_gla = std(log10(nanmean(plt_gla.pow_sens,1)),[],3)/sqrt(size(plt_gla.pow_sens,3));
std_hh  = std(log10(nanmean(plt_hh.pow_sens,1)),[],3)/sqrt(size(plt_hh.pow_sens,3));
std_mue = nanstd(log10(nanmean(plt_mue.pow_sens,1)),[],3)/sqrt(size(plt_mue.pow_sens,3));
std_cnt = nanstd(log10(nanmean(plt_hh_cnt.pow_sens,1)),[],3)/sqrt(size(plt_hh_cnt.pow_sens,3));

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

subplot(4,3,5)
shadedErrorBar(log10(freqoi),par_cnt,std_cnt,{'color',colors(2,:)})
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
par_cnt = mean(nanmean(plt_hh_cnt.corr_sens,1),3);

std_gla = std(nanmean(plt_gla.corr_sens,1),[],3)/sqrt(size(plt_gla.corr_sens,3));
std_hh  = std(nanmean(plt_hh.corr_sens,1),[],3)/sqrt(size(plt_hh.corr_sens,3));
std_mue = std(nanmean(plt_mue.corr_sens,1),[],3)/sqrt(size(plt_mue.corr_sens,3));
std_cnt = std(nanmean(plt_hh_cnt.corr_sens,1),[],3)/sqrt(size(plt_hh_cnt.corr_sens,3));

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

subplot(4,3,5)
shadedErrorBar(log10(freqoi),par_cnt,std_cnt,{'color',colors(2,:)})
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

d_cnt_rest = mean(nanmean(plt_hh_cnt.corr_sens-plt_hh.corr_sens,1),3);
d_cnt_rest_std = std(nanmean(plt_hh_cnt.corr_sens-plt_hh.corr_sens,1),[],3)/sqrt(size(plt_hh_cnt.corr_sens-plt_hh.corr_sens,3));

subplot(4,3,6)
shadedErrorBar(log10(freqoi),d_cnt_rest,d_cnt_rest_std,{'color',colors(2,:)})
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


tmp = squeeze(mean(cat(3,plt_gla.corr_src,plt_hh.corr_src,plt_mue.corr_src),1));
[c,p]=permutest(tmp,zeros(size(squeeze(tmp))),1,0.01,1000,2); h=[c{p<0.05}];

subplot(4,3,4); hold on; box on
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
plot(log10(freqoi(h)),par_all(h),'r.','markersize',8)
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))


par_cnt = mean(nanmean(plt_hh_cnt.corr_src,1),3);
std_cnt = std(nanmean(plt_hh_cnt.corr_src,1),[],3)/sqrt(size(plt_hh_cnt.corr_src,3));
[c,p]=permutest(squeeze(nanmean(plt_hh_cnt.corr_src,1)),zeros(size(squeeze(nanmean(plt_hh_cnt.corr_src,1)))),1,0.01,1000,2); h=[c{p<0.05}];

subplot(4,3,5); hold on; box on
shadedErrorBar(log10(freqoi),par_cnt,std_cnt,{'color',colors(2,:)})
plot(log10(freqoi(h)),par_cnt(h),'k.','markersize',8)
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

subplot(4,3,4); hold on
tmp = squeeze(mean(cat(3,plt_gla.corr_src_df,plt_hh.corr_src_df,plt_mue.corr_src_df),1));
[c,p]=permutest(tmp,zeros(size(squeeze(tmp))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
plot(log10(freqoi(h)),par_all(h),'k.','markersize',8)
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
% title('Derivative: POOLED')

par_cnt = mean(nanmean(plt_hh_cnt.corr_src_df,1),3);
std_cnt = std(nanmean(plt_hh_cnt.corr_src_df,1),[],3)/sqrt(size(plt_hh_cnt.corr_src_df,3));
[c,p]=permutest(squeeze(nanmean(plt_hh_cnt.corr_src_df,1)),zeros(size(squeeze(nanmean(plt_hh_cnt.corr_src_df,1)))),1,0.01,1000,2); h=[c{p<0.05}];

subplot(4,3,5); hold on; box on
shadedErrorBar(log10(freqoi),par_cnt,std_cnt,{'color',colors(2,:)})
plot(log10(freqoi(h)),par_cnt(h),'k.','markersize',8)
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))


print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_derivative_lineplots_v%d.pdf',v))

%% PLOT SOURCE SPACE
% idx_sorted = 1 : 246;
[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');
%
% idx_R = idx_sorted(contains(BNA.tissuelabel(idx_sorted),'Right'));
% idx_L = idx_sorted(contains(BNA.tissuelabel(idx_sorted),'Left'));
% %
% idx_sorted = [idx_R idx_L];

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

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,1);
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

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_correlations_BNA_v%d.pdf',v))

figure_w;
subplot(2,4,5)
[h,p]=ttest(plt_hh_cnt.corr_src_BNA(idx_sorted,:,:),zeros(size(plt_hh_cnt.corr_src_BNA(idx_sorted,:,:))),'dim',3);h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(plt_hh_cnt.corr_src_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,6)

par = plt_hh_cnt.corr_src_BNA(idx_sorted,:,:)-plt_hh.corr_src_BNA(idx_sorted,:,:);
[h,p]=ttest(par,zeros(size(par)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(par,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,7)

[h,p]=ttest(plt_hh_cnt.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_hh_cnt.corr_src_df_BNA(idx_sorted,:,:))),'dim',3); 
imagesc(nanmean(plt_hh_cnt.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,4,8)

[h,p]=ttest(plt_hh_cnt.corr_src_df_BNA(idx_sorted,:,:),plt_hh.corr_src_df_BNA(idx_sorted,:,:),'dim',3); 
imagesc((nanmean(plt_hh_cnt.corr_src_df_BNA(idx_sorted,:,:),3)-nanmean(plt_hh.corr_src_df_BNA(idx_sorted,:,:),3)).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_task_src_correlations_BNA_v%d.pdf',v))

%% PLOT REGIONS OF INTEREST, SENSORY AREAS
load ~/pp/proc/pp_atlas_BNA.mat

is_dt = 0; is_task = 0;

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
subplot(4,3,1); hold on; box off
if ~is_dt && ~is_task
  sig_A1_lr = mean(plt_all.corr_src([A1_l A1_r],:,:),1);
elseif is_dt && ~is_task
  sig_A1_lr = mean(plt_all.corr_src_df([A1_l A1_r],:,:),1);
elseif ~is_dt && is_task
  sig_A1_lr = mean(plt_hh_cnt.corr_src_df([A1_l A1_r],:,:),1);
elseif is_dt && is_task
  sig_A1_lr = mean(plt_hh_cnt.corr_src_df([A1_l A1_r],:,:),1);
end

[c,p]=permutest(squeeze(sig_A1_lr),zeros(size(squeeze(sig_A1_lr))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

subplot(4,3,2); hold on; box off
if ~is_dt && ~is_task
  sig_V1_lr = mean(plt_all.corr_src([V1_l V1_r],:,:),1);
elseif is_dt && ~is_task
  sig_V1_lr = mean(plt_all.corr_src_df([V1_l V1_r],:,:),1);
elseif ~is_dt && is_task
  sig_V1_lr = mean(plt_hh_cnt.corr_src_df([V1_l V1_r],:,:),1);
elseif is_dt && is_task
  sig_V1_lr = mean(plt_hh_cnt.corr_src_df([V1_l V1_r],:,:),1);
end

[c,p]=permutest(squeeze(sig_V1_lr),zeros(size(squeeze(sig_V1_lr))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,3); hold on; box off
if ~is_dt && ~is_task
  sig_M1_lr = mean(plt_all.corr_src([M1_l M1_r],:,:),1);
elseif is_dt && ~is_task
  sig_M1_lr = mean(plt_all.corr_src_df([M1_l M1_r],:,:),1);  
elseif ~is_dt && is_task
  sig_M1_lr = mean(plt_hh_cnt.corr_src_df([M1_l M1_r],:,:),1);
elseif is_dt && is_task
  sig_M1_lr = mean(plt_hh_cnt.corr_src_df([M1_l M1_r],:,:),1);
end

[c,p]=permutest(squeeze(sig_M1_lr),zeros(size(squeeze(sig_M1_lr))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]);% title(sprintf('SOM: r = %.3f',r_M1_lr))

% figure_w;
dlpfc_l=BNA.tissue_5mm==3;
dlpfc_r=BNA.tissue_5mm==4;

subplot(4,3,7); hold on; box off
if ~is_dt && ~is_task
  sig_dlpfc_lr = mean(plt_all.corr_src(find([dlpfc_l+dlpfc_r]),:,:),1);
elseif is_dt && ~is_task
  sig_dlpfc_lr = mean(plt_all.corr_src_df(find([dlpfc_l+dlpfc_r]),:,:),1);
elseif ~is_dt && is_task
  sig_dlpfc_lr = mean(plt_hh_cnt.corr_src_df(find([dlpfc_l+dlpfc_r]),:,:),1);
elseif is_dt && is_task
  sig_dlpfc_lr = mean(plt_hh_cnt.corr_src_df(find([dlpfc_l+dlpfc_r]),:,:),1);
end
[~,p]=ttest(sig_dlpfc_lr,zeros(size(sig_dlpfc_lr)),'dim',3); h=(p*25)<0.05;
shadedErrorBar(log10(freqoi),nanmean(sig_dlpfc_lr,3),std(sig_dlpfc_lr,[],3)/sqrt(size(sig_dlpfc_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_dlpfc_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]);% title(sprintf('V1: r = %.3f',r_V1_lr))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_task%d_dt%d_v%d.pdf',is_task,is_dt,v))


%% PLOT SOURCE MAPS: PUPIL
% load /home/gnolte/meth/templates/mri.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
is_dt = 0; is_task = 1;

for ifoi = [5 11 14 22]
  
  figure_w
  
  if is_dt && ~is_task
    [h,p] = ttest(plt_all.corr_src_df(:,ifoi,:),zeros(size(plt_all.corr_src_df(:,ifoi,:))),'dim',3);
    h=p<(fdr1(p(:),0.05,0));
    par=nanmean(plt_all.corr_src_df(:,ifoi,:),3).*h;
  elseif ~is_dt && ~is_task
    [h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
    h=p<(fdr1(p(:),0.05,0));
    par=nanmean(plt_all.corr_src(:,ifoi,:),3).*h;
  elseif is_dt && is_task
    [h,p] = ttest(plt_hh_cnt.corr_src_df(:,ifoi,:),zeros(size(plt_hh_cnt.corr_src_df(:,ifoi,:))),'dim',3);
    h=p<(fdr1(p(:),0.05,0));
    par=nanmean(plt_hh_cnt.corr_src_df(:,ifoi,:),3).*h;
  elseif ~is_dt && is_task
    [h,p] = ttest(plt_hh_cnt.corr_src(:,ifoi,:),zeros(size(plt_hh_cnt.corr_src(:,ifoi,:))),'dim',3);
    h=p<(fdr1(p(:),0.05,0));
    par=nanmean(plt_hh_cnt.corr_src(:,ifoi,:),3).*h;
  end
%     par(BNA.tissue_5mm>=163 & BNA.tissue_5mm<=174)=0;
%     par(BNA.tissue_5mm>=211 & BNA.tissue_5mm<=214)=0;

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
  
  text(1,1,sprintf('[%.3f %.3f]\n [%.3f Hz]',clim(1),clim(2),ifoi))  

  set(gcf,'renderer','painters')
  print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_avg_task%d_dt%d_f%d_v%d.tiff',is_task,is_dt,ifoi,v))
  
end

%% PLOT CROSS CORRELATION

is_dt = 1;

if is_dt==1
  line_x = 0;
else
  line_x = 0.93;
end

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cols = cbrewer('seq', 'Oranges', 45,'pchip');
cols = cols(1:end-15,:); cols=cols(end:-1:1,:);

figure; set(gcf,'color','w') ;
subplot(2,2,3); hold on; title('Muenster')

for ifreq = 1 : 25
  lags = plt_mue.xcorr_lags{ifreq};
  if is_dt == 0
    plot(lags,nanmean(nanmean(plt_mue.xcorr{ifreq},2),3),'color',cols(ifreq,:))
  else
    plot(lags,nanmean(nanmean(plt_mue.xcorr_df{ifreq},2),3),'color',cols(ifreq,:))
  end
end

axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([line_x line_x],[-0.06 0.04],'color','k')

cols = cbrewer('seq', 'Reds', 35,'pchip');
cols = cols(1:end-5,:); cols=cols(end:-1:1,:);

subplot(2,2,1); hold on; title('Glasgow')

for ifreq = 1 : 25
  lags = plt_gla.xcorr_lags{ifreq};
  if is_dt == 0
    plot(lags,nanmean(nanmean(plt_gla.xcorr{ifreq},2),3),'color',cols(ifreq,:))
  else
    plot(lags,nanmean(nanmean(plt_gla.xcorr_df{ifreq},2),3),'color',cols(ifreq,:))
  end
end

axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([line_x line_x],[-0.06 0.04],'color','k')

cols = cbrewer('seq', 'Blues', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);

subplot(2,2,2); hold on; title('Hamburg')

for ifreq = 1 : 25
  lags = squeeze(plt_hh.xcorr_lags{ifreq});
  if is_dt == 0
    plot(lags,nanmean(nanmean(plt_hh.xcorr{ifreq},2),3),'color',cols(ifreq,:))
  else
    plot(lags,nanmean(nanmean(plt_hh.xcorr_df{ifreq},2),3),'color',cols(ifreq,:))
  end
end
axis([-5 5 -0.06 0.04]); xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots; h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([line_x line_x],[-0.06 0.04],'color','k')

cols = cbrewer('seq', 'Greys', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
for ifreq = 1 : 25
  all.xcorr{ifreq}=cat(2,squeeze(nanmean(plt_gla.xcorr{ifreq},2)),squeeze(nanmean(plt_mue.xcorr{ifreq},2)),squeeze(nanmean(plt_hh.xcorr{ifreq},2)));
  all.xcorr_df{ifreq}=cat(2,squeeze(nanmean(plt_gla.xcorr_df{ifreq},2)),squeeze(nanmean(plt_mue.xcorr_df{ifreq},2)),squeeze(nanmean(plt_hh.xcorr_df{ifreq},2)));
end

subplot(2,2,4); hold on; title('Pooled')

for ifreq = 1 : 25
  if is_dt == 0
    plot(plt_mue.xcorr_lags{ifreq},nanmean(all.xcorr{ifreq},2),'color',cols(ifreq,:))
  else
    plot(plt_mue.xcorr_lags{ifreq},nanmean(all.xcorr_df{ifreq},2),'color',cols(ifreq,:))
  end
end

axis([-5 5 -0.06 0.04])
xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([line_x line_x],[-0.06 0.04],'color','k')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_dt%d_task0_allfreqs_v%d.pdf',is_dt,v))

figure_w
subplot(2,2,4); hold on; title('Pooled')

cols = cbrewer('seq', 'Blues', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);

for ifreq = 1 : 25
  if is_dt == 0
    tmp = nanmean(plt_hh_cnt.xcorr{ifreq},3);
    plot(plt_hh_cnt.xcorr_lags{ifreq},nanmean(tmp,2),'color',cols(ifreq,:))
  else
    tmp = nanmean(plt_hh_cnt.xcorr_df{ifreq},3);
    plot(plt_hh_cnt.xcorr_lags{ifreq},nanmean(tmp,2),'color',cols(ifreq,:))
  end
end

axis([-5 5 -0.06 0.04])
xlabel('Lag [s]'); ylabel('Correlation coeff.');
tp_editplots
h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
line([line_x line_x],[-0.06 0.04],'color','k')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_dt%d_task1_allfreqs_v%d.pdf',is_dt,v))


%% ORDERED XCORR (FRON ANTERIOR TO POSTERIOR)
is_task = 0;
is_dt = 1;

if is_dt==1
  line_x = 0;
else
  line_x = 0.93;
end

cols = cbrewer('seq', 'Greys', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
for ifreq = 1 : 25
  all.xcorr_df_ord{ifreq}=cat(3,plt_gla.xcorr_df_ord{ifreq},plt_mue.xcorr_df_ord{ifreq},plt_hh.xcorr_df_ord{ifreq});
  all.xcorr_ord{ifreq}=cat(3,plt_gla.xcorr_ord{ifreq},plt_mue.xcorr_ord{ifreq},plt_hh.xcorr_ord{ifreq});
end

figure_w;
ff = [1 5; 9 13; 21 25]; color_idx = [1 12 25]; clims = [-0.05 0.05; -0.05 0.05; -0.02 0.02; ]

for iff = 1 : size(ff,1)
  subplot(2,3,iff); hold on
  i = 0; clear sig
  for ifreq = ff(iff,1):ff(iff,2)
    i = i + 1;
    if is_task
      if is_dt==0
        sig(:,:,:,i) = interp1(1:size(plt_hh_cnt.xcorr_ord{ifreq},1),plt_hh_cnt.xcorr_ord{ifreq},linspace(1,size(plt_hh_cnt.xcorr_ord{ifreq},1),size(plt_hh_cnt.xcorr_ord{ff(iff,2)},1)));
      else
        sig(:,:,:,i) = interp1(1:size(plt_hh_cnt.xcorr_df_ord{ifreq},1),plt_hh_cnt.xcorr_df_ord{ifreq},linspace(1,size(plt_hh_cnt.xcorr_df_ord{ifreq},1),size(plt_hh_cnt.xcorr_df_ord{ff(iff,2)},1)));
      end
    else
      if is_dt==0
        sig(:,:,:,i) = interp1(1:size(all.xcorr_ord{ifreq},1),all.xcorr_ord{ifreq},linspace(1,size(all.xcorr_ord{ifreq},1),size(all.xcorr_ord{ff(iff,2)},1)));
      else
        sig(:,:,:,i) = interp1(1:size(all.xcorr_df_ord{ifreq},1),all.xcorr_df_ord{ifreq},linspace(1,size(all.xcorr_df_ord{ifreq},1),size(all.xcorr_df_ord{ff(iff,2)},1)));
      end
    end
  end

  [h,p]=ttest(nanmean(sig,4),zeros(size(nanmean(sig,4))),'dim',3); h = p<fdr1(p(:),0.05,1);
  imagesc(plt_hh.xcorr_lags{ff(iff,2)},1:39,nanmean(nanmean(sig,4),3)'.*h',clims(iff,:))
  set(gca,'ydir','normal','ytick',[],'yticklabels',[]);
  tp_editplots; colormap(cmap)

  axis([-5 5 1 39]);xlabel('Lag [s]'); ylabel('Posterior < Anterior')
  line([line_x line_x],[1 39],'linestyle',':','color',[0.8 0.8 0.8])

end

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_dt%d_task%d_freqsoi_across_sensors_v%d.pdf',is_dt,is_task,v))

%% XCORR FOR SELECTED FREQUENCIES
is_task = 1;
is_dt = 0;

if is_task
  cols = cbrewer('seq', 'Blues', 35,'pchip');
  cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
else
  cols = cbrewer('seq', 'Greys', 35,'pchip');
  cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
end

figure_w

ff = [1 5; 9 13; 21 25]; color_idx = [1 12 25];

for iff = 1 : size(ff,1)
  subplot(2,2,1); hold on

  i = 0; clear sig
  for ifreq = ff(iff,1):ff(iff,2)
    i = i + 1;
    if is_task
      if is_dt
        xcorr_tmp = squeeze(nanmean(plt_hh_cnt.xcorr{ifreq},2));
        sig(:,:,i) = interp1(1:size(xcorr_tmp,1),xcorr_tmp,linspace(1,size(xcorr_tmp,1),size(plt_hh_cnt.xcorr{ff(iff,2)},1)));
      else
        xcorr_tmp = squeeze(nanmean(plt_hh_cnt.xcorr_df{ifreq},2));
        sig(:,:,i) = interp1(1:size(xcorr_tmp,1),xcorr_tmp,linspace(1,size(xcorr_tmp,1),size(plt_hh_cnt.xcorr{ff(iff,2)},1)));
      end
    else
      if is_dt
        sig(:,:,i) = interp1(1:size(all.xcorr_df{ifreq},1),all.xcorr_df{ifreq},linspace(1,size(all.xcorr_df{ifreq},1),size(all.xcorr_df{ff(iff,2)},1)));
      else
        sig(:,:,i) = interp1(1:size(all.xcorr{ifreq},1),all.xcorr{ifreq},linspace(1,size(all.xcorr{ifreq},1),size(all.xcorr{ff(iff,2)},1)));
      end
    end
  end

  sig = zscore(nanmean(nanmean(sig,2),3));
  plot(plt_hh_cnt.xcorr_lags{ff(iff,2)},sig,'color',cols(color_idx(iff),:))
  line([0.93 0.93],[-4 4],'color',[0.8 0.8 0.8],'linestyle',':')
  h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
  axis([-5 5 -6 4]);xlabel('Lag [s]'); tp_editplots; ylabel('Correlation coeff. (z-scored)')
end
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_dt%d_task%d_pooledfreqs_v%d.pdf',is_dt,is_task,v))


%% XCORR: find peak lag and plot spatial distribution

% ifreq = 5;
% % plt_hh.xcorr{ifreq}
% idx = plt_hh.xcorr_lags{ifreq}>=-2 & plt_hh.xcorr_lags{ifreq}<=2;
% lags = plt_hh.xcorr_lags{ifreq}(idx);
% for i =1  : 275
%    for j = 1 : 28
%         [~,max_lag(i,j)] = max(smooth(plt_hh.xcorr{ifreq}(idx,i,j),10));
%         [~,min_lag(i,j)] = min(smooth(nanmean(plt_hh.xcorr{ifreq}(idx,i,j),3),10));
%   end
% end
% para = [];
% para.scale = [-2 2];
% para.cbar = 0;
% figure_w;
% subplot(2,2,1);
% showfield(lags(round(mean(min_lag,2))),[transpose(1:275) lay.pos(1:275,:)],para)
% subplot(2,2,2);
% showfield(lags(round(mean(max_lag,2))),[transpose(1:275) lay.pos(1:275,:)],para)

%% IMGAESC IN HEAD PLOT

task = 0
i = 0; clear sig
for ifreq = 1 : 25
  ifreq
  i = i + 1;
  if task
    sig(:,:,ifreq) = interp1(1:size(hh_cnt.xcorr{ifreq},1),hh_cnt.xcorr{ifreq},linspace(1,size(hh_cnt.xcorr{ifreq},1),size(hh_cnt.xcorr{25},1)));
  else
    if is_dt==0
      sig(:,:,:,ifreq) = interp1(1:size(plt_hh.xcorr{ifreq},1),plt_hh.xcorr{ifreq},linspace(1,size(plt_hh.xcorr{ifreq},1),size(plt_hh.xcorr{25},1)));
    else
      sig(:,:,:,ifreq) = interp1(1:size(plt_hh.xcorr_df{ifreq},1),plt_hh.xcorr_df{ifreq},linspace(1,size(plt_hh.xcorr_df{ifreq},1),size(plt_hh.xcorr_df{25},1)));
    end
  end
end

sig = permute(sig,[2 1 3]);


R=find(startsWith(lay.label(1:275),'MR')); R = R(1:2:end);
L=find(startsWith(lay.label(1:275),'ML')); L = L(1:2:end);
[~,i1]=min(abs(plt_hh.xcorr_lags{25}-(-2)))
[~,i2]=min(abs(plt_hh.xcorr_lags{25}-(2)))
[~,z0] = min(abs(plt_hh.xcorr_lags{25}));

cfg.channel = lay.label([R,L])

sig(:,z0,:)=0;
dummy_freq.powspctrm = sig(:,i1:i2,:);
dummy_freq.time = i1:i2;
% dummy_freq

pars =[];
pars.resolution= 50;
% pars.scale  = [-0.08 0.08];
idx = zeros(275,1); idx(1:4:end)=1;
chan_idx =-1*ones(275,1); chan_idx(2:60:end)=5;
tp_showtfinhead(sig,[chan_idx,lay.pos(1:275,:),10*ones(275,1),ones(275,1)],pars)


%% PLOT SCALING EXPONENTS
v_fooof = 2;
% load fooof results and power spectra
fooof = pp_load_fooof_results(v_fooof);
% fooof.ps_* = power spectra
% fooof.psfit_* = fooof fits of spectra
% fooof.aper_* = offset and slope of 1/f component

% ----------------------------------------
% PLOT OFFSET AND SLOPE OF FOOOF FITS
% ----------------------------------------

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');
%%
figure_w;

data = {'gla';'hh';'mue'};

for idata = 1 : 3
  subplot(4,3,idata); hold on; title(sprintf('Offset - %s',upper(data{idata})))
  eval(sprintf('par_off{idata} = squeeze(nanmean(fooof.aper_%s(1,idx_sorted,3,:)-fooof.aper_%s(1,idx_sorted,1,:),3));',data{idata},data{idata}))
  plot(nanmean(par_off{idata},2),'color',colors(idata,:));
  [h,p] = ttest(par_off{idata},zeros(size(par_off{idata})),'dim',2); h = p<fdr1(p(:),0.05,1);
  if sum(h)>0
    plot(find(h),nanmean(par_off{idata}(find(h),:),2),'r.','markersize',5)
  end
  axis([1 246 -0.1 0.1]); set(gca,'ytick',[-0.1 0 0.1],'yticklabel',[0.1 0 0.1]); 
  line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
  tp_editplots; ylabel('Diff. in offs.'); xlabel('Anterior > Posterior')
  set(gca,'xtick',[],'xticklabel',[]); 

  subplot(4,3,3+idata); hold on; title(sprintf('Slope - %s',upper(data{idata})))
  eval(sprintf('par_slp{idata} = squeeze(nanmean(fooof.aper_%s(2,idx_sorted,3,:)-fooof.aper_%s(2,idx_sorted,1,:),3))',data{idata},data{idata}));
  plot(nanmean(par_slp{idata},2),'color',colors(idata,:));
  [h,p] = ttest(par_slp{idata},zeros(size(par_slp{idata})),'dim',2); h = p<fdr1(p(:),0.05,1);
  if sum(h)>0
    plot(find(h),nanmean(par_slp{idata}(find(h),:),2),'r.','markersize',5)
  end
  axis([1 246 -0.075 0.075]); set(gca,'ytick',[-0.075 0 0.075],'yticklabel',[0.075 0 0.075]); 
  line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
  tp_editplots; ylabel('Diff. in slope'); xlabel('Anterior > Posterior')
  set(gca,'xtick',[],'xticklabel',[]); 

end

subplot(4,3,7); hold on; title(sprintf('Offset - Pooled'))
par_off_all = cat(2,par_off{1},par_off{2},par_off{3});
plot(nanmean(par_off_all,2),'color','k');
[h,p] = ttest(par_off_all,zeros(size(par_off_all)),'dim',2); h = p<fdr1(p(:),0.05,1);
if sum(h)>0
  plot(find(h),nanmean(par_off_all(find(h),:),2),'r.','markersize',5)
end
axis([1 246 -0.05 0.05]); set(gca,'ytick',[-0.05 0 0.05],'yticklabel',[0.05 0 0.05]); 
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in offs.'); xlabel('Anterior > Posterior')
set(gca,'xtick',[],'xticklabel',[]);

subplot(4,3,10); hold on; title(sprintf('Slope - Pooled'))
par_slp_all = cat(2,par_slp{1},par_slp{2},par_slp{3});
plot(nanmean(par_slp_all,2),'color','k');
[h,p] = ttest(par_slp_all,zeros(size(par_slp_all)),'dim',2); h = p<fdr1(p(:),0.05,1);
if sum(h)>0
  plot(find(h),nanmean(par_slp_all(find(h),:),2),'r.','markersize',5)
end
axis([1 246 -0.04 0.04]); set(gca,'ytick',[-0.04 0 0.04],'yticklabel',[0.04 0 0.04]); 
line([0 246],[0 0],'linestyle',':','color',[0.8 0.8 0.8])
tp_editplots; ylabel('Diff. in offs.'); xlabel('Anterior > Posterior')
set(gca,'xtick',[],'xticklabel',[]);

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_fooof_offset_and_slope_v%d.pdf',v_fooof))
%% SPATIAL MAPS OF SLOPE AND OFFSET CHANGE
% START WITH SLOPE
% ---------------------------
[h,p] = ttest(par_slp_all,zeros(size(par_slp_all)),'dim',2); h = p<fdr1(p(:),0.05,1);
idx = find(h);
par_src = zeros(8799,80);
for i = 1:sum(h)
  par_src(BNA.tissue_5mm == idx(i),:) = repmat(par_slp_all(idx(i),:),[sum(BNA.tissue_5mm == idx(i)), 1]);
end

par = nanmean(par_src,2);
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
para = [];
para.colorlimits = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para.colormaps{1} = cmap;
para.orientation = 'axial';
para.dslice_shown = 0.95;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par]); f=get(gcf)
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
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_fooof_slope_spatialmap_v%d.tiff',v_fooof))
% ---------------------------
% OFFSET SECOND
% ---------------------------

[h,p] = ttest(par_off_all,zeros(size(par_off_all)),'dim',2); h = p<fdr1(p(:),0.05,1);
idx = find(h);
par_src = zeros(8799,80);
for i = 1:sum(h)
  par_src(BNA.tissue_5mm == idx(i),:) = repmat(par_off_all(idx(i),:),[sum(BNA.tissue_5mm == idx(i)), 1]);
end

par = nanmean(par_src,2);
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
para = [];
para.colorlimits = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para.colormaps{1} = cmap;
para.orientation = 'axial';
para.dslice_shown = 0.95;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par]); f=get(gcf)
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

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_fooof_offset_spatialmap_v%d.tiff',v_fooof))

%% CORRELATIONS
ffx = 2:0.5:128;
pow = cat(3,fooof.ps_gla,fooof.ps_hh,fooof.ps_mue);
slp  = cat(4,fooof.aper_gla,fooof.aper_hh,fooof.aper_mue);

d_pow = pow(:,:,:,3)-pow(:,:,:,1);
d_slp = squeeze(slp(2,:,3,:)-slp(2,:,1,:));
d_off = squeeze(slp(1,:,3,:)-slp(1,:,1,:));

idx = ~isnan(d_pow(:,1,1));
for isubj = 1 : 28
  tmp1=d_pow(idx,:,isubj);
  tmp2=d_slp(:,isubj);
  tmp3=d_off(:,isubj);
  
  r_powslp(:,isubj)=corr(tmp1',tmp2);
  r_powoff(:,isubj)=corr(tmp1',tmp3);
 
end

[~, p]=ttest(r_powslp,zeros(size(r)),'dim',2); h_slp=p<fdr1(p(:),0.05,1);
[~, p]=ttest(r_powoff,zeros(size(r)),'dim',2); h_off=p<fdr1(p(:),0.05,1);

% ffx = ffx(idx);
figure_w;
par = nan(1,length(ffx)); par(idx) = nanmean(r_powoff,2); 
h = logical(zeros(1,length(ffx))); h(idx) = h_off;
subplot(3,2,1); hold on
plot(log10(ffx),par,'k'); 
plot(log10(ffx(h)),par(:,h),'r.','markersize',8)
axis([log10(2) log10(128) -0.15 0.15]); tp_editplots; box on
xlabel('Frequency [Hz]'); ylabel('Correlation coeff.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

par = nan(1,length(ffx)); par(idx) = nanmean(r_powslp,2); 
h = logical(zeros(1,length(ffx))); h(idx) = h_off;
subplot(3,2,2); hold on
plot(log10(ffx),par,'k'); 
plot(log10(ffx(h)),par(:,h),'r.','markersize',8)
axis([log10(2) log10(128) -0.15 0.15]); tp_editplots; box on
xlabel('Frequency [Hz]'); ylabel('Correlation coeff.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_fooof_correlation_with_changes_in_power_v%d.pdf',v_fooof))