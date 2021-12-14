%%

clear
restoredefaultpath
load ~/pp/proc/pp_atlas_BNA.mat

addpath ~/pconn/matlab/
addpath ~/pp/matlab/
load /home/gnolte/meth/templates/mri.mat
load /home/gnolte/meth/templates/sa_template.mat
addpath ~/Documents/MATLAB/fieldtrip-20190224/
ft_defaults


addpath /home/gnolte/meth/highlevel/

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
v = 2

[plt_gla,plt_hh,plt_mue,plt_all]=pp_load_results(v);

colors = cbrewer('qual', 'Set3', 10,'pchip');
colors = colors(4:6,:);

%% CORRELATION BETWEEN QUADRATIC PARAMETER AND CORR COEFF. - SOURCE SPAC
% REVISION % REVISION % REVISION % REVISION % REVISION % REVISION
% ------------------------------------------

% load results from polynomial fit
load ~/pp/proc/src/pp_invU.mat
freqoi    = 2.^(1:(1/4):7);
pooled_nai = cat(3,zscore(plt_gla.nai_src,1),zscore(plt_hh.nai_src,1),zscore(plt_mue.nai_src,1));

idx = 1:246;

par = zeros(8799,25,81);
for ifreq = 1 : 25
  ifreq
  for isubj = 1 : 81
    for i = 1:246
      par(BNA.tissue_5mm == idx(i),ifreq,isubj) = repmat(p_low22(ifreq,idx(i),1,isubj),[sum(BNA.tissue_5mm == idx(i)), 1]);
    end
  end
end

clear rrr
for ifreq=1:25
  for isubj = 1 : 81
    
    [rrr(ifreq,isubj)]=corr(squeeze(par(:,ifreq,isubj)),squeeze(pooled_nai(:,ifreq,isubj)));
  end
end

rrr = rrr';
[h,p] = ttest(rrr',zeros(size(rrr')),'dim',2); h=p<fdr1(p(:),0.1,0);
s = std(rrr,[],1)/sqrt(size(rrr,1));
figure_w; hold on
shadedErrorBar([],mean(rrr,1),s);
line([0 25],[0 0],'color','k','linestyle',':')
for ii = 1 : max(bwlabel(h)) 
  plot(find(bwlabel(h)==ii),0.1*ones(sum(bwlabel(h)==ii),1),'color',[0.8 0.8 0.8],'linewidth',4);
end
axis([0 25 -0.2 0.1])
tp_editplots
set(gca,'xtick',1:4:25,'xticklabel',freqoi(1:4:end))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_with_invU_v%d.pdf',v))
clear rrr

%% CORRELATION BETWEEN POWER AND CORR COEFF. - SOURCE SPACE
% REVISION % REVISION % REVISION % REVISION % REVISION % REVISION
% ------------------------------------------

freqoi    = 2.^(1:(1/4):7);
pooled_nai = cat(3,zscore(plt_gla.nai_src,1),zscore(plt_hh.nai_src,1),zscore(plt_mue.nai_src,1));
pooled_corr = cat(3,plt_gla.corr_src,plt_hh.corr_src,plt_mue.corr_src);

for isubj=1:81
  for ifreq = 1 : 25
    
    rrr(isubj,ifreq)=corr(pooled_nai(:,ifreq,isubj),pooled_corr(:,ifreq,isubj));
    
  end
end

[h,p] = ttest(rrr',zeros(size(rrr')),'dim',2); h=p<fdr1(p(:),0.1,1);

s = std(rrr,[],1)/sqrt(size(rrr,1));
figure_w; hold on
shadedErrorBar([],mean(rrr,1),s);
line([0 25],[0 0],'color','k','linestyle',':')
for ii = 1 : max(bwlabel(h)) 
  plot(find(bwlabel(h)==ii),0.2*ones(sum(bwlabel(h)==ii),1),'color',[0.8 0.8 0.8],'linewidth',4);
end
tp_editplots
set(gca,'xtick',1:4:25,'xticklabel',freqoi(1:4:end))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_with_power_v%d.pdf',v))

%% REGRESS OUT POWER MAP FROM CORRELATIONS
% REVISION % REVISION % REVISION % REVISION % REVISION % REVISION
% ------------------------------------------

for ifreq = 1 : 25
  ifreq
  for isubj = 1 : 81
  [~,~,pooled_corr_res(:,ifreq,isubj)]=regress(squeeze(pooled_corr(:,ifreq,isubj)),squeeze(pooled_nai(:,ifreq,isubj)));
  end
end

% PLOT IN SOURCE SPACE
% ---------------------------

ifoi = 22;

[h,p]=ttest(pooled_corr_res(:,ifoi,:),zeros(size(pooled_corr_res(:,ifoi,:))),'dim',3); h = p<fdr1(p(:),0.1,0);
par = mean(pooled_corr_res(:,ifoi,:),3).*h;
par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);

clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];

para = [];
para.colorlimits = clim;
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
  
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_residual_task0_dt0_f%d_v%d.tiff',ifoi,v))


%% PLOT POWER SOURCE MAPS
% REVISION % REVISION % REVISION % REVISION % REVISION % REVISION
% ------------------------------------------

ifoi = 22;

% [h,p]=ttest(pooled_nai(:,ifoi,:),zeros(size(pooled_nai(:,ifoi,:))),'dim',3); h = p<fdr1(p(:),0.05,1);
par = mean(pooled_nai(:,ifoi,:),3);
par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);

clim = [min(par) max(par)];

para = [];
para.colorlimits = clim;
para.colormaps{1} = plasma;
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
  
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_power_task0_dt0_f%d_v%d.tiff',ifoi,v))
% -------------------------------

%% PLOT SENSOR LEVEL POWER
% See Supplementary Figure S1

% -----------------
% SENSOR SPACE - POWER
% -----------------
freqoi=2.^(1:(1/4):7);

tot_size = size(plt_mue.corr_sens,3)+size(plt_gla.corr_sens,3)+size(plt_hh.corr_sens,3);

par_gla = mean(log10(nanmean(plt_gla.pow_sens,1)),3);
par_hh  = mean(log10(nanmean(plt_hh.pow_sens,1)),3);
par_mue = mean(log10(nanmean(plt_mue.pow_sens,1)),3);

% par_cnt = nanmean(log10(nanmean(plt_hh_cnt.pow_sens,1)),3);

std_gla = std(log10(nanmean(plt_gla.pow_sens,1)),[],3)/sqrt(size(plt_gla.pow_sens,3));
std_hh  = std(log10(nanmean(plt_hh.pow_sens,1)),[],3)/sqrt(size(plt_hh.pow_sens,3));
std_mue = nanstd(log10(nanmean(plt_mue.pow_sens,1)),[],3)/sqrt(size(plt_mue.pow_sens,3));
% std_cnt = nanstd(log10(nanmean(plt_hh_cnt.pow_sens,1)),[],3)/sqrt(size(plt_hh_cnt.pow_sens,3));

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

all_pow = cat(2,log10(squeeze(nanmean(plt_gla.pow_sens,1))),log10(squeeze(nanmean(plt_hh.pow_sens,1))),log10(squeeze(nanmean(plt_mue.pow_sens,1))))

subplot(4,3,[5 6 8 9])
imagesc(corr(all_pow,all_pow),[0.7 1])
tp_editplots; axis square; colormap(plasma)
set(gca,'xtick',[22 50],'xticklabel',{'HH';'MUE'})
set(gca,'ytick',[22 50],'yticklabel',{'HH';'MUE'})

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_pow_lineplots_v%d.pdf',v))

%% PLOT SOURCE SPACE
% Figure 4B (pupil; v2) and Supplementary Figure S5B (pupil derivative; v1)

load ~/pp/proc/pp_atlas_BNA.mat
% idx_sorted = 1 : 246;
cmap = redblue;
[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

if v==2 
  is_dt = 0;
else
  is_dt = 1;
end

figure_w

if v == 2 

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

  [~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.1,0);
  imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
  tp_editplots; colormap(cmap)
  set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
  xlabel('Frequency [Hz]'); ylabel('Brain region')
 
elseif v==1
  
  subplot(2,4,1)
  [h,p]=ttest(plt_gla.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_gla.corr_src_BNA)),'dim',3);
  imagesc(nanmean(plt_gla.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
  tp_editplots; colormap(cmap)
  set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
  xlabel('Frequency [Hz]'); ylabel('Brain region')

  subplot(2,4,2)
  [h,p]=ttest(plt_hh.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_hh.corr_src_BNA)),'dim',3);
  imagesc(nanmean(plt_hh.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
  tp_editplots; colormap(cmap)
  set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
  xlabel('Frequency [Hz]'); ylabel('Brain region')

  subplot(2,4,3)
  [h,p]=ttest(plt_mue.corr_src_df_BNA(idx_sorted,:,:),zeros(size(plt_mue.corr_src_BNA)),'dim',3);
  imagesc(nanmean(plt_mue.corr_src_df_BNA(idx_sorted,:,:),3).*h,[-0.05 0.05])
  tp_editplots; colormap(cmap)
  set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
  xlabel('Frequency [Hz]'); ylabel('Brain region')

  subplot(2,4,4)

  pooled= cat(3,plt_hh.corr_src_df_BNA(idx_sorted,:,:),plt_gla.corr_src_df_BNA(idx_sorted,:,:),plt_mue.corr_src_df_BNA(idx_sorted,:,:));

  [~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.1,0);
  imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
  tp_editplots; colormap(cmap)
  set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
  xlabel('Frequency [Hz]'); ylabel('Brain region')
else
  error('Version 3? Then change code')
end

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_correlations_BNA_dt%d_v%d.pdf',is_dt,v))

%% PLOT REGIONS OF INTEREST, SENSORY AREAS
% Figure 4A,C (pupil; v=2) and Supplementary Figure S5A,C (pupil derivative; v==1)
% -----------------------------

% load atlas
load ~/pp/proc/pp_atlas_BNA.mat

if v == 2 || v==3 
  is_dt = 0; 
  is_task = 0; % task is always 0
else
  is_dt = 1; 
  is_task = 0;
end

V1_r = find(BNA.tissue_5mm==205);
V1_l = find(BNA.tissue_5mm==206);
A1_r = find(BNA.tissue_5mm==71);
A1_l = find(BNA.tissue_5mm==72);
M1_r = find(BNA.tissue_5mm==159);
M1_l = find(BNA.tissue_5mm==160);
ACC_l= find(BNA.tissue_5mm==179);
ACC_r= find(BNA.tissue_5mm==180);

figure_w;

subplot(4,3,1); hold on; box off
if ~is_dt && ~is_task
  sig_pooled = mean(plt_all.corr_src(:,:,:),1);
elseif is_dt && ~is_task
  sig_pooled = mean(plt_all.corr_src_df(:,:,:),1);
end

[c,p]=ttest(squeeze(sig_pooled),zeros(size(squeeze(sig_pooled))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_pooled),zeros(size(squeeze(sig_pooled))),1,0.01,1000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_pooled,3),std(sig_pooled,[],3)/sqrt(size(sig_pooled,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_pooled(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,2); hold on; box off
if ~is_dt && ~is_task
  sig_V1_lr = mean(plt_all.corr_src([V1_l V1_r],:,:),1);
elseif is_dt && ~is_task
  sig_V1_lr = mean(plt_all.corr_src_df([V1_l V1_r],:,:),1);
end

[c,p]=ttest(squeeze(sig_V1_lr),zeros(size(squeeze(sig_V1_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_V1_lr),zeros(size(squeeze(sig_V1_lr))),1,0.01,1000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,5); hold on; box off
if ~is_dt && ~is_task
  sig_A1_lr = mean(plt_all.corr_src([A1_l A1_r],:,:),1);
elseif is_dt && ~is_task
  sig_A1_lr = mean(plt_all.corr_src_df([A1_l A1_r],:,:),1);
end

[c,p]=ttest(squeeze(sig_A1_lr),zeros(size(squeeze(sig_A1_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_A1_lr),zeros(size(squeeze(sig_A1_lr))),1,0.01,1000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
end

line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

subplot(4,3,8); hold on; box off
if ~is_dt && ~is_task
  sig_M1_lr = mean(plt_all.corr_src([M1_l M1_r],:,:),1);
elseif is_dt && ~is_task
  sig_M1_lr = mean(plt_all.corr_src_df([M1_l M1_r],:,:),1);  
end

[c,p]=ttest(squeeze(sig_M1_lr),zeros(size(squeeze(sig_M1_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_M1_lr),zeros(size(squeeze(sig_M1_lr))),1,0.01,1000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]);% title(sprintf('SOM: r = %.3f',r_M1_lr))


subplot(4,3,11); hold on; box off
if ~is_dt && ~is_task
  sig_dlpfc_lr = mean(plt_all.corr_src([ACC_l ACC_r],:,:),1);
elseif is_dt && ~is_task
  sig_dlpfc_lr = mean(plt_all.corr_src_df(find([ACC_l ACC_r]),:,:),1);
end
% 
[c,p]=ttest(squeeze(sig_dlpfc_lr),zeros(size(squeeze(sig_dlpfc_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_dlpfc_lr),zeros(size(squeeze(sig_dlpfc_lr))),1,0.01,1000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_dlpfc_lr,3),std(sig_dlpfc_lr,[],3)/sqrt(size(sig_dlpfc_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_dlpfc_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]);% title(sprintf('V1: r = %.3f',r_V1_lr))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_task%d_dt%d_v%d.pdf',is_task,is_dt,v))


%% PLOT SOURCE MAPS: PUPIL
% Figure 4D (pupil; v=2) and Supplementary Figure 5D (pupil derivative; v=3)
% --------------------------------------------------

% pooled / gla (glasgow) / hh (hamburg)  / mue (muenster)
site_to_plot = 'gla';

% load /home/gnolte/meth/templates/mri.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

% cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
if v==2 || v==3 
 is_dt = 0; 
 is_task = 0;
else
 is_dt = 1; 
 is_task = 0; 
end

cmap = redblue;


for ifoi = [5 11 14 22]
  
  figure_w
  
  if strcmp(site_to_plot,'pooled')
    if is_dt && ~is_task
      [h,p] = ttest(plt_all.corr_src_df(:,ifoi,:),zeros(size(plt_all.corr_src_df(:,ifoi,:))),'dim',3);
      p_adj=(fdr1(p(:),0.1,0));
      h=p<(fdr1(p(:),0.1,0));
      par=nanmean(plt_all.corr_src_df(:,ifoi,:),3).*h;
    elseif ~is_dt && ~is_task
      [h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
      h=p<(fdr1(p(:),0.1,0));
      p_adj=(fdr1(p(:),0.1,0));
      par=nanmean(plt_all.corr_src(:,ifoi,:),3).*h;
    end  
  elseif strcmp(site_to_plot,'gla') % Plot Glasgow data
    if is_dt && ~is_task
      par=nanmean(plt_gla.corr_src_df(:,ifoi,:),3);
    elseif ~is_dt && ~is_task
      par=nanmean(plt_gla.corr_src(:,ifoi,:),3);
    end  
  elseif strcmp(site_to_plot,'hh') % Plot Hamburg data
    if is_dt && ~is_task
      par=nanmean(plt_hh.corr_src_df(:,ifoi,:),3);
    elseif ~is_dt && ~is_task
      par=nanmean(plt_hh.corr_src(:,ifoi,:),3);
    end  
  elseif strcmp(site_to_plot,'mue') % Plot Muenster data
    if is_dt && ~is_task
      par=nanmean(plt_mue.corr_src_df(:,ifoi,:),3);
    elseif ~is_dt && ~is_task
      par=nanmean(plt_mue.corr_src(:,ifoi,:),3);
    end  
  end
  % project onto fine grid
  par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);
  
  mask=h;
  save(sprintf('~/pp/proc/src/pp_sourcemaps_mask_f%d_v%d.mat',ifoi,v),'mask')

  
  clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
%   clim = clim-0.1*clim;
  para = [];
  para.colorlimits = clim;
  para.colormaps{1} = cmap;
  para.orientation = 'axial';
  if strcmp(site_to_plot,'pooled') 
    para.thresh = 1
  else
    para.thresh = 0;
  end
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
  
  text(1,1,sprintf('[%.3f %.3f]\n [%.3f Hz] | p_adj = %.4f',clim(1),clim(2),ifoi,p_adj))  

  set(gcf,'renderer','painters')
  
  if strcmp(site_to_plot,'pooled')
    print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_avg_task%d_dt%d_f%d_v%d.tiff',is_task,is_dt,ifoi,v))
  else
    print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_avg_%s_task%d_dt%d_f%d_v%d.tiff',site_to_plot,is_task,is_dt,ifoi,v))    
  end
end

%% PLOT CROSS CORRELATION
% !!!!!!!IMPORTANT!!!!!!
% Cross correlations should be plotted for v==1, as this is the no-lag
% analysis. In v2, pupils are already shifted with respect to MEG
% !!!!!!!!!

if v == 2 || v == 3
  error('Makes no sense for v2/v3')
end

is_dt = 0;

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

subplot(2,2,2); hold on; title('Hamburg')

for ifreq = 1 : 25
  lags = squeeze(plt_hh_cnt.xcorr_lags{ifreq});
  if is_dt == 0
    plot(lags,nanmean(nanmean(plt_hh_cnt.xcorr{ifreq},2),3),'color',cols(ifreq,:))
  else
    plot(lags,nanmean(nanmean(plt_hh_cnt.xcorr_df{ifreq},2),3),'color',cols(ifreq,:))
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


%% ORDERED XCORR (FRON ANTERIOR TO POSTERIOR)
is_task = 0;
is_dt = 1;

to_plot = 'all';

if is_dt==1
  line_x = 0;
else
  line_x = 0.93;
end

cols = cbrewer('seq', 'Greys', 35,'pchip');
cols = cols(3:end-3,:); cols=cols(end:-1:1,:);
for ifreq = 1 : 25
  if strcmp(to_plot,'all')
    all.xcorr_df_ord{ifreq}=cat(3,plt_gla.xcorr_df_ord{ifreq},plt_mue.xcorr_df_ord{ifreq},plt_hh.xcorr_df_ord{ifreq});
    all.xcorr_ord{ifreq}=cat(3,plt_gla.xcorr_ord{ifreq},plt_mue.xcorr_ord{ifreq},plt_hh.xcorr_ord{ifreq});
  elseif strcmp(to_plot,'gla')
    all.xcorr_df_ord{ifreq}=plt_gla.xcorr_df_ord{ifreq};
    all.xcorr_ord{ifreq}=plt_gla.xcorr_ord{ifreq};
  elseif strcmp(to_plot,'hh')
    all.xcorr_df_ord{ifreq}=plt_hh.xcorr_df_ord{ifreq};
    all.xcorr_ord{ifreq}=plt_hh.xcorr_ord{ifreq};
  elseif strcmp(to_plot,'mue')
    all.xcorr_df_ord{ifreq}=plt_mue.xcorr_df_ord{ifreq};
    all.xcorr_ord{ifreq}=plt_mue.xcorr_ord{ifreq};
  end
end
figure_w;
ff = [1 5; 9 13; 21 25]; color_idx = [1 12 25]; clims = [-0.05 0.05; -0.05 0.05; -0.02 0.02; ]

for iff = 1 : size(ff,1)
  subplot(2,3,iff); hold on
  i = 0; clear sig
  for ifreq = ff(iff,1):ff(iff,2)
    i = i + 1;
    
    if is_dt==0
      sig(:,:,:,i) = interp1(1:size(all.xcorr_ord{ifreq},1),all.xcorr_ord{ifreq},linspace(1,size(all.xcorr_ord{ifreq},1),size(all.xcorr_ord{ff(iff,2)},1)));
    else
      sig(:,:,:,i) = interp1(1:size(all.xcorr_df_ord{ifreq},1),all.xcorr_df_ord{ifreq},linspace(1,size(all.xcorr_df_ord{ifreq},1),size(all.xcorr_df_ord{ff(iff,2)},1)));
    end
  end
  
  [h,p]=ttest(nanmean(sig,4),zeros(size(nanmean(sig,4))),'dim',3); 
  if strcmp(to_plot,'all')
    h = p<fdr1(p(:),0.1,0);
  end
  imagesc(plt_hh.xcorr_lags{ff(iff,2)},1:39,nanmean(nanmean(sig,4),3)'.*h',clims(iff,:))
  set(gca,'ydir','normal','ytick',[],'yticklabels',[]);
  tp_editplots; colormap(cmap)
  
  axis([-5 5 1 39]);xlabel('Lag [s]'); ylabel('Posterior < Anterior')
  line([line_x line_x],[1 39],'linestyle',':','color',[0.8 0.8 0.8])
  
end

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_%s_xcorr_dt%d_task%d_freqsoi_across_sensors_v%d.pdf',to_plot,is_dt,is_task,v))

%% XCORR FOR SELECTED FREQUENCIES
is_task = 0;
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
  
  sig_stats = nanmean(sig,3);
  [h,p]=ttest(sig_stats,zeros(size(sig_stats)),'dim',2);
  
%   real_lags = plt_hh_cnt.xcorr_lags{ff(iff,2)}(plt_hh_cnt.xcorr_lags{ff(iff,2)}>-1.5 & plt_hh_cnt.xcorr_lags{ff(iff,2)} < 1.5);
%   [i,j]=min(sig(plt_hh_cnt.xcorr_lags{ff(iff,2)}>-1.5 & plt_hh_cnt.xcorr_lags{ff(iff,2)} < 1.5,:)); 
%   all_lags(iff,:)=real_lags(j);
  
  sig = zscore(nanmean(nanmean(sig,2),3));
  plot(plt_hh_cnt.xcorr_lags{ff(iff,2)},sig,'color',cols(color_idx(iff),:))
  line([0.93 0.93],[-6 4],'color',[0 0 0],'linestyle',':')
  line([0 0],[-6 4],'color',[0.8 0.8 0.8],'linestyle',':')
  h=colorbar; colormap(gca,cols); h.Label.String = 'Frequency [Hz]'; h.TickLabels={2;128}; h.Ticks=[0 1];
  axis([-5 5 -6 4]);xlabel('Lag [s]'); tp_editplots; ylabel('Correlation coeff. (z-scored)')
end
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_xcorr_dt%d_task%d_pooledfreqs_v%d.pdf',is_dt,is_task,v))







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


%% COMPARISON OF PUPIL AND PUPIL DERIVATIVE: PLOT SOURCE SPACE
% --------------------------------------------

load ~/pp/proc/pp_atlas_BNA.mat
% idx_sorted = 1 : 246;
cmap = redblue;
[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

figure_w

subplot(2,4,4)

pooled = plt_all1.corr_src_df_BNA-plt_all2.corr_src_BNA;

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.1,0);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
xlabel('Frequency [Hz]'); ylabel('Brain region')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_correlations_BNA_dt1_vs_dt0_v%d.pdf',v))

%% PLOT REGIONS OF INTEREST, SENSORY AREAS
load ~/pp/proc/pp_atlas_BNA.mat

M1 = [-42 -26 54; 38 -32 48];
V1 = [-20 -86 18; 16 -80 26];
A1 = [-54 -22 10; 52 -24 12];

V1_r = find(BNA.tissue_5mm==205);
V1_l = find(BNA.tissue_5mm==206);
A1_r = find(BNA.tissue_5mm==71);
A1_l = find(BNA.tissue_5mm==72);
M1_r = find(BNA.tissue_5mm==159);
M1_l = find(BNA.tissue_5mm==160);
ACC_l= find(BNA.tissue_5mm==179);
ACC_r= find(BNA.tissue_5mm==180);

figure_w;

subplot(4,3,1); hold on; box off
sig_pooled = mean(plt_all1.corr_src_df(:,:,:),1)-mean(plt_all2.corr_src(:,:,:),1);

[c,p]=ttest(squeeze(sig_pooled),zeros(size(squeeze(sig_pooled))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_pooled),zeros(size(squeeze(sig_pooled))),1,0.01,10000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_pooled,3),std(sig_pooled,[],3)/sqrt(size(sig_pooled,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_pooled(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,2); hold on; box off
sig_V1_lr = mean(plt_all1.corr_src_df([V1_l V1_r],:,:),1)-mean(plt_all2.corr_src([V1_l V1_r],:,:),1);


[c,p]=ttest(squeeze(sig_V1_lr),zeros(size(squeeze(sig_V1_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_V1_lr),zeros(size(squeeze(sig_V1_lr))),1,0.01,10000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_V1_lr,3),std(sig_V1_lr,[],3)/sqrt(size(sig_V1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_V1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('V1: r = %.3f',r_V1_lr))

subplot(4,3,5); hold on; box off

sig_A1_lr = mean(plt_all1.corr_src_df([A1_l A1_r],:,:),1)-mean(plt_all2.corr_src([A1_l A1_r],:,:),1);

[c,p]=ttest(squeeze(sig_A1_lr),zeros(size(squeeze(sig_A1_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_A1_lr),zeros(size(squeeze(sig_A1_lr))),1,0.05,10000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_A1_lr,3),std(sig_A1_lr,[],3)/sqrt(size(sig_A1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_A1_lr(:,h,:),3),'k.','markersize',8)
end

line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

subplot(4,3,8); hold on; box off

sig_M1_lr = mean(plt_all1.corr_src_df([M1_l M1_r],:,:),1)-mean(plt_all2.corr_src([M1_l M1_r],:,:),1);

[c,p]=ttest(squeeze(sig_M1_lr),zeros(size(squeeze(sig_M1_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_M1_lr),zeros(size(squeeze(sig_M1_lr))),1,0.01,10000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_M1_lr,3),std(sig_M1_lr,[],3)/sqrt(size(sig_M1_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_M1_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]);% title(sprintf('SOM: r = %.3f',r_M1_lr))


subplot(4,3,11); hold on; box off

sig_dlpfc_lr = mean(plt_all1.corr_src_df([ACC_l ACC_r],:,:),1)-mean(plt_all2.corr_src([ACC_l ACC_r],:,:),1);

[c,p]=ttest(squeeze(sig_dlpfc_lr),zeros(size(squeeze(sig_dlpfc_lr))),'dim',2); h=p<fdr1(p(:),0.1,0);

% [c,p]=permutest(squeeze(sig_dlpfc_lr),zeros(size(squeeze(sig_dlpfc_lr))),1,0.01,10000,1); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(sig_dlpfc_lr,3),std(sig_dlpfc_lr,[],3)/sqrt(size(sig_dlpfc_lr,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(sig_dlpfc_lr(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Correlation')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]);% title(sprintf('V1: r = %.3f',r_V1_lr))

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_dt1_vs_dt0_v%d.pdf',v))


%% PLOT SOURCE MAPS: PUPIL
% load /home/gnolte/meth/templates/mri.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = redblue;

for ifoi = [5 11 14 22]
  
  figure_w

  pooled = plt_all1.corr_src_df(:,ifoi,:) -plt_all2.corr_src(:,ifoi,:);
  
    [h,p] = ttest(pooled,zeros(size(pooled)),'dim',3);
    
    h=p<(fdr1(p(:),0.1,0));
    par=nanmean(pooled,3).*h;
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
  print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_avg_dt1_vs_dt0_f%d_v%d.tiff',ifoi,v))
  
end

%% PLOT INVERTED U


% Load data computed in powerspectra.m
fooof22 = pp_load_fooof_results(2);
fooof11 = pp_load_fooof_results(1);

tmp22 = cat(4,nanmean(fooof22.pxx_seg_hh,5),fooof22.pxx_seg_gla,fooof22.pxx_seg_mue);
tmp11 = cat(4,nanmean(fooof11.pxx_seg_dt_hh,5),fooof11.pxx_seg_dt_gla,fooof11.pxx_seg_dt_mue);


fxx = fooof22.fxx;
fooof22.fxx = fxx;

% set notch filter range to nan
idx = (fooof22.fxx>=47 & fooof22.fxx<=53) |  (fooof22.fxx>=97 & fooof22.fxx<=103);
tmp22(idx,:,:,:) = nan;
tmp11(idx,:,:,:) = nan;

pooled22 = zeros(size(tmp22),'double');
pooled11 = zeros(size(tmp11),'double');

% interpolate nans (speed up at some point)
for ireg = 1 : size(pooled22,2)
  ireg
  for iseg = 1 : size(pooled22,3)
    for isubj = 1 : size(pooled22,4)
      
      pooled22(:,ireg,iseg,isubj) = fillmissing(tmp22(:,ireg,iseg,isubj),'spline');
      pooled11(:,ireg,iseg,isubj) = fillmissing(tmp11(:,ireg,iseg,isubj),'spline');
      
    end
  end
end

%%

clear mean_dat22 mean_dat11 pooled_f11

pooled_f22 = zeros(25,size(pooled22,2),size(pooled22,3),size(pooled22,4),'double');
pooled_f11 = zeros(25,size(pooled11,2),size(pooled11,3),size(pooled11,4),'double');

f = 2.^(1:0.25:7);
for ii = 1 : 25
  [wavelet, outp]= tp_mkwavelet(f(ii),0.5,400,0);
  freq_range(ii,:) = outp.freq_range;
  pooled_f22(ii,:,:,:) = squeeze(mean(pooled22(fxx>=freq_range(ii,1) & fxx<=freq_range(ii,2),:,:,:),1));
  pooled_f11(ii,:,:,:) = squeeze(mean(pooled11(fxx>freq_range(ii,1) & fxx<freq_range(ii,2),:,:,:),1));
end

x = 1 :14;

p_low22 = zeros(size(pooled_f22,1),size(pooled_f22,2),3,size(pooled_f22,4),'double');
p_low11 = zeros(size(pooled_f22,1),size(pooled_f22,2),3,size(pooled_f22,4),'double');

for isubj = 1  : size(pooled_f22,4)
  
  %   m_subj(isubj) = mean(mean(mean1/(pooled_f22(:,:,:,isubj),1),2),3);
  isubj
  rand_idx = randperm(20);
  for iff = 1  : 25
    for i = 1 : 246
      
      dat = squeeze(pooled_f22(iff,i,:,isubj))';
      
      % POLYFIT WITH DEG 2 and
      if any(isnan(dat))
        sum(isnan(dat))
        dat = zscore(fillmissing(dat,'spline'));
      else
        dat = zscore(dat);
      end
      
      p_low22(iff,i,:,isubj)=polyfit(x,dat,2);
      mean_dat22(iff,i,:,isubj) = dat;
      
      % PERMUTE PUPIL SEGMENTS
      dat = squeeze(pooled_f11(iff,i,:,isubj))';
      % POLYFIT WITH DEG 2 and
      if any(isnan(dat))
        % interpolation is only needed if pupil segments are not ordered by
        % quantiles
        dat = zscore(fillmissing(dat,'spline'));
      else
        dat = zscore(dat);
      end
      
      p_low11(iff,i,:,isubj)=polyfit(x,dat,2);
      
      mean_dat11(iff,i,:,isubj) = dat;
      
      clear dat
      
    end
  end
end

save(['~/pp/proc/src/pp_invU.mat'],'p_low11','p_low22')
%%
close

v=2;
conjunction = 0;
for ifoi = [11 14]
  
  h = ones(246,1);
  idx = find(h);
  %
  if v==2
    ppp = p_low22(:,:,1,:);%./mean(pooled_f(:,:,:,:),3);
  else
    ppp = p_low11(:,:,1,:);%./mean(pooled_f(:,:,:,:),3);
  end
  %
  % load(sprintf('~/pvals_%d.mat',ifoi),'h')
  % mm = h; clear h
  [h,p]=ttest(ppp(ifoi,:,:,:),zeros(size(ppp(ifoi,:,:,:))),'dim',4);
  h = p<fdr1(p(:),0.1,0);
  
  masked = 1;
  % prc_max = mean(abs(10.5-idx_max),3);
  if masked == 1
    prc_max = mean(ppp(ifoi,:,:,:),4).*h;
  else
    prc_max = mean(ppp(ifoi,:,:,:),4);
  end
  
  
  par = zeros(8799,1);
  hh=zeros(8799,1);
  % mask=zeros(8799,1);
  for i = 1:246
    par(BNA.tissue_5mm == idx(i),:) = repmat(prc_max(:,idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
    hh(BNA.tissue_5mm == idx(i),:) = repmat(h(idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
    %   mask(BNA.tissue_5mm == idx(i),:) = repmat(mm(idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
  end
  
  if conjunction==1
    par = zeros(8799,1);
    %   par(mask&(~hh))=30; iii = 5;
    % par(mask&hh)=30; iii = 1;
    par((~mask)&hh)=30; iii = 3;
    
  end
  sum(par==1)
  if conjunction==1
    cmap = cbrewer('qual', 'Set1', 5,'pchip');
    cmap = cmap(iii,:);
    %   cmap = cmap;
  else
    cmap = redblue;
  end
  para = [];
  para.colorlimits = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
  % para.colorlimits = [min(par(:)) max(par(:))-5];
  % para.colorlimits = [5.1 8]
  tmp=plasma;
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
  
  text(1,1,sprintf('[%.8f %.3f]\n [%.8f Hz]',para.colorlimits(1),para.colorlimits(2),ifoi))
  
  set(gcf,'renderer','painters')
  if ~conjunction
    print(gcf,'-dpdf',sprintf('~/pp/plots/pp_quad_masked%d_spatialmap_f%d_v%d.tiff',masked,ifoi,v))
  else
    print(gcf,'-dpdf',sprintf('~/pp/plots/pp_quad_spatialmap_conj%d_f%d_v%d.tiff',iii,ifoi,v))
  end
end

%%

ifoi = 22;

[h,p]=ttest(squeeze(p_low22(ifoi,:,1,:)),zeros(size(squeeze(p_low22(ifoi,:,1,:)))),'dim',2);
h=p<fdr1(p,0.1,0);
% h = p<0.01;
par = squeeze(mean(p_low22(ifoi,:,1,:),4));
idx = 1:246;
par1 = zeros(8799,1);
for i = 1:246
  par1(BNA.tissue_5mm == idx(i)) = repmat(par(:,idx(i)).*h(idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
end

% [h,p]=ttest(pooled,zeros(size(pooled)),'dim',3);
% par = mean(pooled(:,ifoi,:),3).*(p(:,ifoi)<0.0005);

cmap      = redblue;
para      = [];
para.clim = [-max([abs([min(par1(:)) max(par1(:))])]) max([abs([min(par1(:)) max(par1(:))])])];
para.cmap = cmap;
para.grid = BNA.grid_5mm/10;
para.dd   = 0.5;
para.fn   = sprintf('~/test_ACC.png');
tp_plot_surface(par1,para)

%%

v = 1;
freqs = 2.^(1:0.25:7);

figure_w

load ~/pp/proc/pp_atlas_BNA.mat

if v == 2
  is_dt = 0;
  is_task = 0;
else
  is_dt = 1;
  is_task = 0;
end


figure_w


subplot(4,3,1); hold on; box off

if v==2
  ppp = p_low22(:,:,1,:);
  sig_V1 = squeeze(mean(p_low22(:,[205 206],1,:),2));
  sig_A1 = squeeze(mean(p_low22(:,[71 72],1,:),2));
  sig_M1 = squeeze(mean(p_low22(:,[159 160],1,:),2));
  sig_ACC = squeeze(mean(p_low22(:,[179 180],1,:),2));
elseif v==1
  ppp = p_low11(:,:,1,:);
  sig_V1 = squeeze(mean(p_low11(:,[205 206],1,:),2));
  sig_A1 = squeeze(mean(p_low11(:,[71 72],1,:),2));
  sig_M1 = squeeze(mean(p_low11(:,[159 160],1,:),2));
  sig_ACC = squeeze(mean(p_low11(:,[179 180],1,:),2));
end

par = squeeze(mean(ppp,2));
[h,p]=ttest(par,zeros(size(par)),'dim',2)
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(par,2),std(par,[],2)/sqrt(size(par,2)))
plot(log10(freqs(find(p<0.05))),mean(par(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(par(find(p<0.01),:),2),'w.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(par(find(p<p_fdr),:),2),'y.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
axis([.3 2.11 -0.013 0.013])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))


subplot(4,3,2); hold on; box off
[h,p]=ttest(sig_V1,zeros(size(sig_V1)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_V1,2),std(par,[],2)/sqrt(size(sig_V1,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_V1(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_V1(find(p<0.01),:),2),'w.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_V1(find(p<p_fdr),:),2),'y.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])

subplot(4,3,5); hold on; box off
[h,p]=ttest(sig_A1,zeros(size(sig_A1)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_A1,2),std(par,[],2)/sqrt(size(sig_A1,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_A1(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_A1(find(p<0.01),:),2),'w.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_A1(find(p<p_fdr),:),2),'y.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])

subplot(4,3,8); hold on; box off
[h,p]=ttest(sig_M1,zeros(size(sig_M1)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_M1,2),std(par,[],2)/sqrt(size(sig_M1,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_M1(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_M1(find(p<0.01),:),2),'w.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_M1(find(p<p_fdr),:),2),'y.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])

subplot(4,3,11);hold on; box off
[h,p]=ttest(sig_ACC,zeros(size(sig_ACC)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_ACC,2),std(par,[],2)/sqrt(size(sig_ACC,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_ACC(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_ACC(find(p<0.01),:),2),'w.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_ACC(find(p<p_fdr),:),2),'y.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])
% close al
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_invertedU_dt%d_v%d.pdf',is_dt,v))


%%
clear m22  poly_par22 m11 poly_par11
% is_dt=0;
% pooled11=double(cat(4,double(nanmean(fooof11.pxx_seg_dt_hh,5)),double(fooof11.pxx_seg_dt_gla),double(fooof11.pxx_seg_dt_mue)));

% tmp22 = squeeze(mean(mean(pooled22,4),2));
% tmp11 = squeeze(mean(mean(pooled11,4),2));

for i=1:81
  i
  for ii=1:246
    fxx = 2:0.5:128;
    idx_low=freqs>=2 & freqs<=4;
    idx_alpha=freqs>=8 & freqs<=16;
    idx_high=freqs>=64 & freqs<=128;
    
    m22(1,ii,i,:)=mean(mean_dat22(idx_low,ii,:,i),1);
    m22(2,ii,i,:)=mean(mean_dat22(idx_alpha,ii,:,i),1);
    m22(3,ii,i,:)=mean(mean_dat22(idx_high,ii,:,i),1);
    
    m11(1,ii,i,:)=mean(mean_dat11(idx_low,ii,:,i),1);
    m11(2,ii,i,:)=mean(mean_dat11(idx_alpha,ii,:,i),1);
    m11(3,ii,i,:)=mean(mean_dat11(idx_high,ii,:,i),1);
    
  end
end
%
poly_par22(1,:) = polyfit(1:size(m22,4),squeeze(mean(mean(m22(1,:,:,:),2),3)),2);
poly_par22(2,:) = polyfit(1:size(m22,4),squeeze(mean(mean(m22(2,:,:,:),2),3)),2);
poly_par22(3,:) = polyfit(1:size(m22,4),squeeze(mean(mean(m22(3,:,:,:),2),3)),2);

poly_par11(1,:) = polyfit(1:size(m11,4),squeeze(mean(mean(m11(1,:,:,:),2),3)),2);
poly_par11(2,:) = polyfit(1:size(m11,4),squeeze(mean(mean(m11(2,:,:,:),2),3)),2);
poly_par11(3,:) = polyfit(1:size(m11,4),squeeze(mean(mean(m11(3,:,:,:),2),3)),2);

figure_w
subplot(3,3,1); hold on
plot(squeeze(mean(mean(m22(1,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:size(m22,4),poly_par22(1,1).*(1:size(m22,4)).^2 + poly_par22(1,2).*(1:size(m22,4)) + poly_par22(1,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,2); hold on
plot(squeeze(mean(mean(m22(2,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:size(m22,4),poly_par22(2,1).*(1:size(m22,4)).^2 + poly_par22(2,2).*(1:size(m22,4)) + poly_par22(2,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,3); hold on
plot(squeeze(mean(mean(m22(3,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:size(m22,4),poly_par22(3,1).*(1:size(m22,4)).^2 + poly_par22(3,2).*(1:size(m22,4)) + poly_par22(3,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.95 1.03])

%
subplot(3,3,4); hold on
plot(squeeze(mean(mean(m11(1,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:25,poly_par11(1,1).*(1:25).^2 + poly_par11(1,2).*(1:25) + poly_par11(1,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,5); hold on
plot(squeeze(mean(mean(m11(2,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:25,poly_par11(2,1).*(1:25).^2 + poly_par11(2,2).*(1:25) + poly_par11(2,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,6); hold on
plot(squeeze(mean(mean(m11(3,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:25,poly_par11(3,1).*(1:25).^2 + poly_par11(3,2).*(1:25) + poly_par11(3,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_invU_avg_all.pdf'))

%% PLOT SPECTRAL EXPONENT AND FOOOF ANALYSIS


fooof22 = pp_load_fooof_results(2);
fooof11 = pp_load_fooof_results(1);
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%%
load ~/pp/proc/pp_atlas_BNA.mat
[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');
colors = cbrewer('qual', 'Set3', 10,'pchip');
colors = colors(4:6,:);



%% SLOPES PER SEGMENT
 
pooled22 = squeeze(cat(5,nanmean(fooof22.slopes_seg_hh,6),fooof22.slopes_seg_gla,fooof22.slopes_seg_mue));
pooled11 = squeeze(cat(5,nanmean(fooof11.slopes_seg_hh_dt,6),fooof11.slopes_seg_gla_dt,fooof11.slopes_seg_mue_dt));

r22 = corr((1:14)',mean(nanmean(pooled22,3),1)');
r11 = corr((1:14)',mean(nanmean(pooled11,3),1)');

figure_w;

subplot(3,4,1); hold on
plot(1:14,squeeze(mean(nanmean(pooled22,3),1)),'.','markersize',20)
lsline
tp_editplots; xlabel('Pupil bin'); ylabel('Spectral exponent')
axis([-1 15 0.90 0.98])
text(0,0.91,sprintf('r=%.4f',r22),'fontsize',6)

% GLASGOW
% --------
subplot(3,4,2); hold on
plot(1:14,squeeze(mean(nanmean(fooof22.slopes_seg_gla,5),1)),'.','markersize',20,'color',colors(1,:))
plot(1:14,squeeze(mean(nanmean(nanmean(fooof22.slopes_seg_hh,6),5),1)),'.','markersize',20,'color',colors(2,:))
plot(1:14,squeeze(mean(nanmean(fooof22.slopes_seg_mue,5),1)),'r.','markersize',20,'color',colors(3,:))

lsline
tp_editplots; xlabel('Pupil bin'); ylabel('Spectral exponent')
axis([-1 15 0.82 1.02])


subplot(3,4,5); hold on

plot(1:14,squeeze(mean(nanmean(pooled11,3),1)),'.','markersize',20)
lsline

tp_editplots; xlabel('Pupil bin'); ylabel('Spectral exponent')
axis([-1 15 0.92 1])
text(0,0.925,sprintf('r=%.4f',r11),'fontsize',6)

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_slope_vs_pupil.pdf'))

%%
v=1
close

V1_r = find(BNA.tissue_5mm==205);
V1_l = find(BNA.tissue_5mm==206);
A1_r = find(BNA.tissue_5mm==71);
A1_l = find(BNA.tissue_5mm==72);
M1_r = find(BNA.tissue_5mm==159);
M1_l = find(BNA.tissue_5mm==160);
ACC_l= find(BNA.tissue_5mm==179);
ACC_r= find(BNA.tissue_5mm==180);

cmap = cbrewer('seq', 'YlGnBu', 10,'pchip'); 
cmap = cmap(end:-1:1,:);

if v==2
pooled = squeeze(cat(2,fooof22.slope_hh,fooof22.slope_gla,fooof22.slope_mue));
elseif v==1
pooled = squeeze(cat(2,fooof11.slope_df_hh,fooof11.slope_df_gla,fooof11.slope_df_mue));  
end

figure_w
subplot(5,3,[1 2 4 5]); hold on

r = (rand(81,1)-0.5)/4;

plot(ones(81,1)+r,squeeze(mean(pooled(:,:),1)),'.','markersize',12,'color',[0.7 0.7 0.7])
% plot(1.4,squeeze(mean(mean(pooled(:,:),1),2)),'k.','markersize',20)
[~,p]=ttest(squeeze(mean(pooled(:,:),1)))

plot(2*ones(81,1)+r,squeeze(mean(pooled([205, 206],:),1)),'.','markersize',12,'color',cmap(5,:))
% plot(2.4,squeeze(mean(mean(pooled([205, 206],:),1),2)),'.','markersize',20)

[~,p1]=ttest(squeeze(mean(pooled([205, 206],:),1)))

plot(3*ones(81,1)+r,squeeze(mean(pooled([71, 72],:),1)),'.','markersize',12,'color',cmap(6,:))
[~,p2]=ttest(squeeze(mean(pooled([71, 72],:),1)))

plot(4*ones(81,1)+r,squeeze(mean(pooled([159, 160],:),1)),'.','markersize',12,'color',cmap(7,:))
[~,p3]=ttest(squeeze(mean(pooled([159, 160],:),1)))

plot(5*ones(81,1)+r,squeeze(mean(pooled([179, 180],:),1)),'.','markersize',12,'color',cmap(8,:))
[~,p4]=ttest(squeeze(mean(pooled([179, 180],:),1)))

axis([0.5 5.5 -0.40 0.40])

line([0.8 1.2],[mean(mean(pooled(:,:),1),2) mean(mean(pooled(:,:),1),2)],'linewidth',2)
line([1.8 2.2],[mean(mean(pooled([205, 206],:),1),2) mean(mean(pooled([205, 206],:),1),2)],'linewidth',2)
line([2.8 3.2],[mean(mean(pooled([71, 72],:),1),2) mean(mean(pooled([71, 72],:),1),2)],'linewidth',2)
line([3.8 4.2],[mean(mean(pooled([159, 160],:),1),2) mean(mean(pooled([159, 160],:),1),2)],'linewidth',2)
line([4.8 5.2],[mean(mean(pooled([179, 180],:),1),2) mean(mean(pooled([179, 180],:),1),2)],'linewidth',2)

line([0.5 4.5],[0 0],'linewidth',0.5,'color',[0.7 0.7 0.7],'linestyle','--')

tp_editplots

set(gca,'xtick',1:5,'xticklabels',{'Average';'Visual cortex';'Auditory cortex';'Somatosensory cortex';'Anterior cigulate cortex'})
set(gca,'ytick',[-0.8 -0.4 0 0.4 0.8],'yticklabels',[-0.8 -0.4 0 0.4 0.8])
ylabel('Correlation')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_slope_ROIs_v%d.pdf',v))

%%

addpath /home/gnolte/meth/highlevel/
load /home/gnolte/meth/templates/sa_template;

load /home/gnolte/meth/templates/mri.mat
figure_w
is_slope = 1; is_dt = 1;


if is_slope && ~is_dt
  v = 2;
  pooled = cat(2,fooof22.slope_gla,fooof22.slope_hh,fooof22.slope_mue);
elseif ~is_slope && ~is_dt
  v = 2;
  pooled = cat(2,fooof22.offset_gla,fooof22.offset_hh,fooof22.offset_mue);
elseif is_slope && is_dt
  v = 1;
  pooled = cat(2,fooof11.slope_df_gla,fooof11.slope_df_hh,fooof11.slope_df_mue);
elseif ~is_slope && is_dt
  v = 1;
  pooled = cat(2,fooof11.offset_df_gla,fooof11.offset_df_hh,fooof11.offset_df_mue);
end

cmap = redblue;
par = zeros(8799,81);

for i = 1 : size(pooled,1)
  par(BNA.tissue_5mm==i,:)=repmat(pooled(i,:),[sum(BNA.tissue_5mm==i) 1]); 
end

[h,p]=ttest(par,zeros(size(par)),'dim',2); h = p < fdr1(p(:),0.1,0);

par = nanmean(par,2).*h;
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

text(1,1,sprintf('[%.3f %.3f]',clim(1),clim(2)))

set(gcf,'renderer','painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_fooof_sourcemap_dt%d_slp%d_v%d.tiff',is_dt,is_slope,v))


%% EXAMPLE FOOOF
close
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
isubj = 16;
iblock = 1;
v = 1;

load(sprintf('~/pp/proc/src/pp_gla_collected_fooof_s%d_b%d_v%d.mat',isubj,iblock,v))
load(sprintf('~/pp/proc/src/pp_gla_src_powerspectra_s%d_b%d_v%d.mat',isubj,iblock,v))
 
pup = pup(~isnan(pxx(1,1,:)));
pxx=pxx(:,:,~isnan(pxx(1,1,:)));

[r,i]=min(corr(squeeze(aper(2,:,:))',pup(~isnan(pup))'));

pxx = pxx(:,:,~isnan(pup));
% 
[idx1] = pup>prctile(pup,75);
[idx2] = pup<prctile(pup,25);

% idx1 = 3;
% idx2 = 100;

ap1 = aper(:,i,idx1);
aper1 = (-ap1(2,:)'.*log10(3:0.5:40)+repmat(ap1(1),[75 1])')';
aper1 = mean(aper1,2);

gg1 = squeeze(mean(gg(i,:,idx1),3));

ap2 = aper(:,i,idx2);
aper2 = (-ap2(2,:)'.*log10(3:0.5:40)+repmat(ap2(1),[75 1])')';
aper2 = mean(aper2,2);

figure_w; hold on
% plot(3:0.5:40,gg1,'k'); hold on;
% plot(3:0.5:40,gg2);
% plot(log10(3:0.5:40),aper1,'k:')
plot(log10(3:0.5:40),mean(mean(pxx(3:77,:,idx1),2),3),'k:')
% plot(log10(3:0.5:40),aper2,'r-')
plot(log10(3:0.5:40),mean(mean(pxx(3:77,:,idx2),2),3),'r-')

%% TAKE fitted 1/f distribution out of empirical spectra
freqoi=2.^(1:(1/4):7); ffx = 2:0.5:128;
colors = cbrewer('qual', 'Set3', 10,'pchip'); colors = colors(4:6,:);

clear ps_hh ps_mue ps_gla ps_hh_df ps_mue_df ps_gla_df

for ifreq = 1 : length(freqoi)
  [wavelet,opt]=tp_mkwavelet(freqoi(ifreq),0.5,400,0.8);
  idx = find(ffx>=opt.freq_range(1) & ffx<=opt.freq_range(2));
  ps_hh(:,ifreq,:) = squeeze(nanmean(nanmean(fooof22.ps_hh_corrected(:,idx,:,:),2),4));
  ps_gla(:,ifreq,:) = squeeze(nanmean(fooof22.ps_gla_corrected(:,idx,:),2));
  ps_mue(:,ifreq,:) = squeeze(nanmean(fooof22.ps_mue_corrected(:,idx,:),2));
  ps_hh_df(:,ifreq,:) = squeeze(nanmean(nanmean(fooof11.ps_hh_df_corrected(:,idx,:,:),2),4));
  ps_gla_df(:,ifreq,:) = squeeze(nanmean(fooof11.ps_gla_df_corrected(:,idx,:),2));
  ps_mue_df(:,ifreq,:) = squeeze(nanmean(fooof11.ps_mue_df_corrected(:,idx,:),2));
end

figure_w;
subplot(4,4,[1]); hold on
par = squeeze(nanmean(ps_gla(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); 
h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(1,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(1,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[2]); hold on
par = squeeze(nanmean(ps_hh(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2)); 
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(2,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(2,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04]);
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[3]); hold on
par = squeeze(nanmean(ps_mue(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(3,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(3,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[4]); hold on
par = squeeze(nanmean(cat(3,ps_gla,ps_hh,ps_mue),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color','k'})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color','k')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[5]); hold on
par = squeeze(nanmean(ps_gla_df(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(1,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(1,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[6]); hold on
par = squeeze(nanmean(ps_hh_df(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(2,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(2,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04]);
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[7]); hold on
par = squeeze(nanmean(ps_mue_df(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(3,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(3,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[8]); hold on
par = squeeze(nanmean(cat(3,ps_gla_df,ps_hh_df,ps_mue_df),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[~,p]=ttest(par,zeros(size(squeeze(par))),'dim',2); h = p<fdr1(p(:),0.1,0);
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color','k'})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color','k')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_figure100d.pdf',v))
