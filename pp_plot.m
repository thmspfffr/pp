%%

clear
restoredefaultpath

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load /home/gnolte/meth/templates/mri.mat

ft_defaults


addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

[plt_gla,plt_hh,plt_all]=pp_load_results();

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

pooled = cat(3,plt_hh.corr_sens_ord,plt_gla.corr_sens_ord);

subplot(2,3,3);
imagesc(nanmean(pooled,3),[-0.03 0.03])
colormap(cmap); tp_editplots; axis square;
set(gca,'ydir','normal','xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_anterior_post_v%d.pdf',v))

%% PLOT AVERAGES ACROSS SPACE, HAMBURG, GLASGOW, ALL

% -----------------
% SENSOR SPACE - POWER
% -----------------
freqoi=2.^(1:(1/4):7); 


par_gla = mean(log10(nanmean(plt_gla.pow_sens,1)),3);
par_hh = mean(log10(plt_hh.pow_sens),2);
std_gla = std(log10(nanmean(plt_gla.pow_sens,1)),[],3)/sqrt(size(plt_gla.pow_sens,3));
std_hh = std(log10(plt_hh.pow_sens),[],2)/sqrt(size(plt_hh.pow_sens,2));

par_all = mean(cat(2,log10(squeeze(nanmean(plt_gla.pow_sens,1))),log10(plt_hh.pow_sens)),2);
std_all = std(cat(2,log10(squeeze(nanmean(plt_gla.pow_sens,1))),log10(plt_hh.pow_sens)),[],2)/sqrt(size(plt_hh.pow_sens,2));

figure_w
subplot(4,3,1)
shadedErrorBar(log10(freqoi),par_gla,std_gla,'b')
axis([.3 2.11 -24 -21])
% line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,2)
shadedErrorBar(log10(freqoi),par_hh,std_hh,'r')
axis([.3 2.11 -24 -21])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,3)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -24 -21])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Log-Power')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))


% -----------------
% SENSOR SPACE - CORRELATIONS
% -----------------
par_gla = mean(nanmean(plt_gla.corr_sens,1),3);
par_hh = mean(plt_hh.corr_sens,2);
std_gla = std(nanmean(plt_gla.corr_sens,1),[],3)/sqrt(size(plt_gla.corr_sens,3));
std_hh = std(plt_hh.corr_sens,[],2)/sqrt(size(plt_hh.corr_sens,2));

par_all = mean(cat(2,squeeze(nanmean(plt_gla.corr_sens,1)),plt_hh.corr_sens),2);
std_all = std(cat(2,squeeze(nanmean(plt_gla.corr_sens,1)),plt_hh.corr_sens),[],2)/sqrt(50);



subplot(4,3,7)
shadedErrorBar(log10(freqoi),par_gla,std_gla,'b')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,8)
shadedErrorBar(log10(freqoi),par_hh,std_hh,'r')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,9)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

% -----------------
% SOURCE SPACE
% -----------------
par_gla = mean(nanmean(plt_gla.corr_src,1),3);
par_hh = mean(mean(plt_hh.corr_src,1),3);
std_gla = std(nanmean(plt_gla.corr_src,1),[],3)/sqrt(size(plt_gla.corr_src,3));
std_hh = std(nanmean(plt_hh.corr_src,1),[],3)/sqrt(size(plt_hh.corr_src,2));

par_all = nanmean(mean(cat(3,plt_gla.corr_src,plt_hh.corr_src),1),3);
std_all = std(nanmean(cat(3,plt_gla.corr_src,plt_hh.corr_src),1),[],3)/sqrt(50);


subplot(4,3,10)
shadedErrorBar(log10(freqoi),par_gla,std_gla,'b')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))


subplot(4,3,11)
shadedErrorBar(log10(freqoi),par_hh,std_hh,'r')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

subplot(4,3,12)
shadedErrorBar(log10(freqoi),par_all,std_all,'k')
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_lineplots_v%d.pdf',v))

%% PLOT SOURCE SPACE - FROM ANTERIOR TO POSTERIOR
figure_w
subplot(2,3,1)
[h,p]=ttest(plt_gla.corr_src_BNA,zeros(size(plt_gla.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_gla.corr_src_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,2)
[h,p]=ttest(plt_hh.corr_src_BNA,zeros(size(plt_hh.corr_src_BNA)),'dim',3); 
imagesc(nanmean(plt_hh.corr_src_BNA,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

subplot(2,3,3)

pooled= cat(3,plt_hh.corr_src_BNA,plt_gla.corr_src_BNA); 
[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,0);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_correlations_BNA_v%d.pdf',v))

%% PLOT REGIONS OF INTEREST, SENSORY AREAS

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

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_v%d.pdf',v))

%% PLOT NAI / POWER
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

for ifoi = 12
% ifoi = 14;

[h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
h=p<(fdr1(p(:),0.05,0));
par=nanmean(plt_all.corr_src(:,ifoi,:),3).*h;

clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.75;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par])
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_avg_f%d_v%d.tiff',ifoi,v))

end

%% PLOT CROSS FREQUENCY




