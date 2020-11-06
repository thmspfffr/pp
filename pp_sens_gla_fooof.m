%% pp_sens_gla_fooof
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 4
% -------------------------
v = 3;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
% -------------------------

addpath('~/Documents/MATLAB/fieldtrip-20181231/')
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

addpath /home/gnolte/meth/highlevel/
%%
% -------------------------
for isubj = 1:24
  
  clear data dat pupil pup 
  
  for iblock = 1:1
    %
    fn = sprintf('pp_sens_gla_fooof_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pupil = data.trial{1}';
      f_sample = data.fsample;
      
      % load meg data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
      end
      
      
      artifPnts=data.cfg.artfctdef.visual.artifact;
      
      cfg=[];
      cfg.layout='4D248.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:248),data.label);
      
    
     
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
%     pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
%     pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    
    pupil_df = diff(pupil);
   
    data.avg = data.trial{1}'; %data.trial{1} = [];
                
    opt.n_win = 4000; % 10s segment length, i.e., 0.1:0.1:100
    opt.n_shift = 4000; % no overlap
    
    nseg=floor((size(data.avg,1)-opt.n_win)/opt.n_shift+1);
    clear pxx fxx pup pup_df
    ff = 2:0.5:40;
     
%     pxx = nan(size(ff,2),248,nseg);
    for iseg = 1 : nseg
        
        seg_dat = data.avg((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win,:);
        seg_pup = mean(pupil((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
        [pxx(:,:,iseg),fxx]=pwelch(seg_dat,hanning(800),[],ff,400,'power');
        pup(iseg) = seg_pup;
        seg_pup = mean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
        pup_df(iseg) = seg_pup;
    end


    save([outdir fn '.mat'],'pxx','fxx','pup','pup_df')
    tp_parallel(fn,outdir,0)
    
    clear outp clear pxx pup_df fxx seg_pup 
    
  end
end

error('!')
%% LOAD EXPONENTS

v=4

r_all = nan(248,24);
for isubj = 1 : 24
    
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d_v%d.mat',isubj,v))
    load(sprintf('~/pp/proc/src/chanidx_s%d.mat',isubj))
    r_all(chanidx>0,isubj) = r(chanidx(chanidx>0));
    
%     outp.tp_sens_pow(outp.,ifreq) = tmp();
    
end




cfg=[];
cfg.layout='4D248.lay';
lay = ft_prepare_layout(cfg);
minmax_gla=[min(lay.pos(:,2)) max(lay.pos(:,2))];
ser_gla = linspace(minmax_gla(1),minmax_gla(2),40);
r_ord= zeros(size(ser_gla,2)-1,24);

for i = 1 : size(ser_gla,2)-1
  idx = lay.pos(:,2)<ser_gla(i+1) & lay.pos(:,2)>ser_gla(i);
  r_ord(i,:,:) = nanmean(r_all(idx,:),1);
end




pars            = [];
pars.markersize = 0;
pars.linewidth  = 9;
pars.cbar       = 0;
pars.scale      = [-0.1 0.1]
pars.cmap       = jet;
pars.resolution = 600;

figure; set (gcf,'color','w')

par = nanmean(r_all,2);
par(isnan(par)) = nanmean(par);
showfield_colormap(par.*h,lay.pos(1:248,:),pars);
drawnow

%% SOURCE SPACE

clear r_all
for isubj = 1 : 24
    isubj
    try
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_exp_s%d_v%d.mat',isubj,v));
%     load(sprintf('/home````r/tpfeffer/pp/proc/src/pp_sens_gla_fooof_s%d_b1_v3.mat',isubj))
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_slp_s%d_v%d.mat',isubj,v));

    load(sprintf('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_s%d_b1_v%d.mat',isubj,v));
    r_all(:,isubj) = r(trans);
    pup_df = pup_df(~isnan(pup_df));
    r_all_df(:,isubj)=corr(slp',pup_df');
    r_all_df(:,isubj) = r_all_df(trans,isubj);
    
    catch me
        r_all(:,isubj) = nan(8799,1);
        r_all_df(:,isubj) = nan(8799,1);
    end
    
%     outp.tp_sens_pow(outp.,ifreq) = tmp();
    
end

load /home/gnolte/meth/templates/mri.mat
load /home/gnolte/meth/templates/sa_template.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% [h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
% h=p<(fdr1(p(:),0.05,0));
par=nanmean(r_all_df(:,:),2)
% par=r;
% par_stats=squeeze(nanmean(nanmean(plt_hh.cf_corr(:,eval(sprintf('meg_f%d',ifreq)),eval(sprintf('pup_f%d',ifreq)),:),3),2));
% [h,p]=ttest(par_stats,zeros(size(par_stats)),'dim',2); h = p<(fdr1(p(:),0.05,0));
% par=par.*h;
% 
clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.75;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par])
set(gcf,'renderer','painters')
% print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_crossfreq_f%d_v%d.tiff',ifreq,v))

%% HAMBURG

v= 4
clear r_all
for isubj = SUBJLIST
    isubj
    try
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_hh_src_fooof_exp_s%d.mat',isubj,v));
%     load(sprintf('/home````r/tpfeffer/pp/proc/src/pp_sens_gla_fooof_s%d_b1_v3.mat',isubj))
    load(sprintf('/home/tpfeffer/pp/proc/src/pp_hh_src_fooof_slp_s%d.mat',isubj,v));

%     load(sprintf('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_s%d_b1_v%d.mat',isubj,v));
    r_all(:,isubj) = r(trans);
    pup_df = pup_df(~isnan(pup_df));
    r_all_df(:,isubj)=corr(slp',pup_df');
    r_all_df(:,isubj) = r_all_df(trans,isubj);
    
    catch me
        r_all(:,isubj) = nan(8799,1);
        r_all_df(:,isubj) = nan(8799,1);
    end
    
%     outp.tp_sens_pow(outp.,ifreq) = tmp();
    
end

load /home/gnolte/meth/templates/mri.mat
load /home/gnolte/meth/templates/sa_template.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% [h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
% h=p<(fdr1(p(:),0.05,0));
par=nanmean(r_all_df(:,:),2)
% par=r;
% par_stats=squeeze(nanmean(nanmean(plt_hh.cf_corr(:,eval(sprintf('meg_f%d',ifreq)),eval(sprintf('pup_f%d',ifreq)),:),3),2));
% [h,p]=ttest(par_stats,zeros(size(par_stats)),'dim',2); h = p<(fdr1(p(:),0.05,0));
% par=par.*h;
% 
clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.75;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par])
set(gcf,'renderer','painters')
% print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_corr_sourcemap_crossfreq_f%d_v%d.tiff',ifreq,v))





