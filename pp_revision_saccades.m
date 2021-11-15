%% pp_revision_saccades
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
v = 1;
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
lag = 0;
% -------------------------
% VERSION 3: with pupil lag
% -------------------------
% v = 2;
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% lag = 1;
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

freqoi=2.^(1:(1/4):7); % 2-128 Hz as per Hipp et al. (2012) Nat Neurosci

%% START WITH HAMBURG DATA
% -------------------------
subj_counter = 0;

for isubj = 25:34
  subj_counter=subj_counter+1;
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
      % data from 'pp_prepare_data.m'
      load(sprintf('~/pp/data/ham/pp_rest_s%d_b%d_v%d.mat',isubj,iblock,1))
      load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
    catch me
      disp('Sth went wrong')
      continue
    end
    
    cfg=[];
    cfg.layout='CTF275.lay';
    lay = ft_prepare_layout(cfg);
    [~, outp.chanidx] = ismember(lay.label(1:275),label(startsWith(label,'M')));
    
    % bp-filter and resample pupil
    % ------
    k = 2; f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005; hil_lo = 2;
    hil_Wn = [hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    tmp = filtfilt(bhil, ahil, pupil(:,[1 4]));
    pupil(:,[1 4])=tmp;
    pupil = resample(pupil,400,1000);
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------
    
    pupil = pupil(end:-1:1,:);
    dat = dat(:,end:-1:1);
    
    len = min([size(pupil,1) size(dat,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    dat = dat(:,1:len);
    pupil = pupil(1:len,:);
    
    dat = dat(:,end:-1:1);
    pupil = pupil(end:-1:1,:);
    %     art_idx = [diff(pupil(:,6)>0.1)==1; zeros(1,1)];
    %     blink_idx = [diff(pupil(:,7)>0.1)==1; zeros(1,1)];
    %     art_idx = [find(diff(art_idx)==1) find(diff(art_idx)==-1)];
    %     blink_idx = [find(diff(blink_idx)==1) find(diff(blink_idx)==-1)];
    % ------
    saccs=tp_detect_microsaccades(pupil(:,2:3),400,5,6);
    
    %     rH(isubj, iblock) = corr(pupil(:,2),pupil(:,4));
    %     rV(isubj, iblock) = corr(pupil(:,3),pupil(:,4));
    %
    %     [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(1,isubj) = k(i);
    %     [n,k]=hist(pupil(:,3),200); [~,i]=max(n); fix(2,isubj) = k(i);
    %
    %     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
    %     dist = sqrt(((pupil(:,2)-fix(1,isubj)).^2) + ((pupil(:,3)-fix(2,isubj)).^2) );
    %
    %     rExc(isubj,iblock) = corr(dist(~isnan(pupil(:,3))),pupil(~isnan(pupil(:,3)),4));
    %
    pupil = zscore(pupil(:,4));
    
    tmp = []; meg = []; k = 0;
    for i = 1 : size(saccs,1)
      if (saccs(i,1)+800)>size(pupil,1)
        tmp(1:1001,i)=nan;
      elseif (saccs(i,1)-200)<1
        tmp(1:1001,i)=nan;
      else
        k = k +1;
        tmp(:,i) = pupil(saccs(i,1)-200:saccs(i,1)+800);
        meg(:,:,k) = dat(:,saccs(i,1)-200:saccs(i,1)+800);
      end
    end
    pupil_locked=nanmean(tmp,2);
    
    % COMPUTE FFT
    % ------------
    segleng = 200;
    segshift = 20;
    nseg=floor((size(meg,2)-segleng)/segshift+1);
    t = -100:20:1001;
    pxx = zeros(101,size(meg,1),nseg,'single');
    tmp_csd = zeros(size(meg,1),size(meg,1),101);
    pupil_count = 0; csd_count = 0; pxx_counter = 0;
    pup_seg = zeros(41,1);
    for isacc = 1 : size(meg,3)
      pupil_count = pupil_count+1; 
      isacc
      tmp_pxx = nan(101,size(meg,1),nseg);
      for iseg = 1 : nseg
        meg_seg  = squeeze(meg(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng,isacc));
        win = (hanning(200)./sum(hanning(200)));
        tmp_pup = squeeze(pupil((iseg-1)*segshift+1:(iseg-1)*segshift+segleng));
        pup_seg(iseg)  = pup_seg(iseg)+sum(win.*tmp_pup);
        
        if any(isnan(squeeze(meg_seg(1,:))))
          continue
        else
          [tmp_pxx(:,:,iseg),fxx] = pwelch(meg_seg',hanning(200),[],200,400,'power'); 
        end
      end
      if any(isnan(tmp_pxx(1,1,:)))
        continue
      else
        pxx_counter=pxx_counter+1;
        pxx = pxx + tmp_pxx;
      end
    end
    
    pup_seg = pup_seg ./ pupil_count; 
    
    for ifreq = 1 : size(pxx,1)
    	csd(:,:,ifreq) = sqrt(squeeze(pxx(ifreq,:,:)))*sqrt(squeeze(pxx(ifreq,:,:)))';
    end
  
%     
    pxx_src = zeros(nseg,8799,101); 
    nai_src = zeros(nseg,8799,101); 
    for ifreq = 1 : 101
      ifreq
      para          = [];
      para.reg      = 0.05;
      filt = tp_beamformer(csd(:,:,ifreq),sa.L_genemaps_aal,para);
      tmp_pxx = sqrt(squeeze(pxx(ifreq,:,:)));
      pxx_src(:,:,ifreq) = (tmp_pxx'*filt).^2;
      
      para          = [];
      para.reg      = 0.05;
      Lr = reshape(sa.L_genemaps_aal,[size(sa.L_genemaps_aal,1) 8799*3]); csd_noise = Lr*Lr';
      [~,noise]      = tp_beamformer(real(csd_noise),sa.L_genemaps_aal,para);
      nai_src(:,:,ifreq) =  pxx_src(:,:,ifreq)./repmat(noise,[1 size(pxx_src,1)])';
    end
    
    save([outdir sprintf('pp_revision_saccades_TFR_hh_isubj%d_iblock%d.mat',isubj,iblock)],'pup_seg','fxx','pxx_src','pupil_locked','nai_src','pxx','-v7.3')

    clear pxx_src csd nai_src  pup_seg csd
    
  end
end

error('!')

 
subj_counter = 0; clear pxx_all pxx_src_all
for isubj = 25 : 34
subj_counter=subj_counter+1
load([outdir sprintf('pp_revision_saccades_TFR_hh_isubj%d_iblock%d.mat',isubj,iblock)])

pxx_all(:,:,subj_counter,iblock)=100*(squeeze(nanmean(pxx(2:end,:,:),2))-squeeze(nanmean(pxx(2:end,:,1),2)))./squeeze(nanmean(pxx(2:end,:,1),2));
pxx_src_all(:,:,subj_counter,iblock)=squeeze(100.*(nanmean(pxx_src(:,:,2:end),2)-nanmean(pxx_src(1,:,2:end),2))./nanmean(pxx_src(1,:,2:end),2));

end

figure_w;

imagesc(nanmean(nanmean(pxx_all,3),4),[-2.5 2.5]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 101],'color','k','linestyle','--')
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:100,'yticklabel',fxx(2:5:end))
colormap(redblue); axis([1 41 1 61])
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_TFR_hh.pdf'))


% mean_pup_seg = squeeze(nanmean(nanmean(pup_seg(:,:,25:34,:),4),1));
% mean_pxx = nanmean(all_pxx(:,:,25:34,:),4);
load /home/gnolte/meth/templates/mri.mat
%%
par = 100*(nanmean(pxx_src_all,4)-nanmean(pxx_src_all(1,:,:,:),4))./nanmean(pxx_src_all(1,:,:,:),4);
par = par(10,:,find(fxx==10));
bl = log10(nanmean(all_pxx(1,:,fxx==10,:),4));
  para = [];
  para.colorlimits = [-30 -27];
  para.colormaps{1} = redblue;
  para.orientation = 'axial';

  para.dslice_shown = 0.95;
  para.colorbar= 0;
  
  tp_showmri_transp(mri,para,[sa.grid_BNA_5mm bl']); f=get(gcf);
  
% all_pxx
  
%%
% for ifreq = 1 : 101
%   for isubj = 1 : 10
%     rr(ifreq,isubj)=corr(mean_pup_seg(:,isubj),squeeze(mean_pxx(ifreq,:,isubj))');
%     
%     
%   end
% end

% plot(mean(rr,2))

  
% figure_w; hold on
% % all_tmp = zscore(nanmean(all_tmp(:,25:34,:),3));
% all_tmp = nanmean(all_tmp(:,25:34,:),3);
% all_tmp = all_tmp-nanmean(all_tmp(1:100,:),1);
% s = std(all_tmp,[],2)/sqrt(size(all_tmp,2));
% shadedErrorBar(1:901,nanmean(all_tmp,2),s)
% line([0 901],[0 0],'color','k','linestyle',':')
% line([100 100],[-2 2],'color','k','linestyle',':')
% xlabel('Time [ms]'); ylabel('Pupil size')
% tp_editplots; axis([0 800 -120 120])
% set(gcf,'color','w')
% set(gca,'xtick',[100 300 500 700 901],'xticklabel',[0 500 1000 1500 2000])
%
%
% print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_hamburg.pdf'))

figure_w;

imagesc(mean_pxx,[-0.025 0.025]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 101],'color','k','linestyle','--')
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:101,'yticklabel',fxx(1:5:end))
colormap(redblue)
% conc_dat = all_tmp;

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_TFR_hh.pdf'))

pup_seg = nanmean(all_tmp(:,25:34,:),3);

pxx = nanmean(nanmean(all_pxx(:,:,25:34,:),3),4);

%% GLASGOW DATA
% -------------------------
% subj2,3,5, 21: impossible
    % subj6,, doable?
rExc = []; clear all_tmp all_pxx n_saccs

subj_counter = 0;
for isubj = 1:24
  if any(isubj==[5,9])
    continue
  else
    subj_counter = subj_counter + 1;
  end
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    %     try
    % load pupil data
    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
    end
%         saccs = data.cfg.artfctdef.sacc.artifact;

    for ii = 2 : size(data.cfg.artfctdef.blink.artifact,1)-1
      data.trial{1}(1:2,data.cfg.artfctdef.blink.artifact(ii,1)-150:data.cfg.artfctdef.blink.artifact(ii,2)+150) = nan;
    end
    
    saccs = data.cfg.artfctdef.sacc.artifact;
    pupil = data.trial{1}';
    f_sample = data.fsample;
    
    % load meg data
    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
    end
    
    artifPnts=data.cfg.artfctdef.visual.artifact;
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    tmp = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil(:,4)=tmp;
        
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------
    
    dat = data.trial{1};  clear data
    
    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
    
    for iart = 1 : size(artifPnts,1)
      dat(:,artifPnts(iart,1):artifPnts(iart,2))=NaN;
    end
    
    clear tp_csd
    pupil(:,1:2)=interp1(1:size(pupil(:,1),1),pupil(:,1:2),1:size(pupil(:,1),1),'pchip');
    saccs=tp_detect_microsaccades(pupil(:,1:2),f_sample,5,6);

%     rH(subj_counter) = corr(pupil(:,1),pupil(:,4));
%     rV(subj_counter) = corr(pupil(:,2),pupil(:,4));
%     
%     [n,k]=hist(pupil(:,1),200); [~,i]=max(n); fix(1,subj_counter) = k(i);
%     [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(2,subj_counter) = k(i);
%     
    %     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
%     dist = sqrt(((pupil(:,1)-fix(1,subj_counter)).^2) + ((pupil(:,2)-fix(2,subj_counter)).^2) );
%     
%     rExc(subj_counter) = corr(dist,pupil(:,4));
    
    pupil = zscore(pupil(:,4));
    
    tmp = []; meg = []; k = 0;
    for i = 1 : size(saccs,1)
      if (saccs(i,1)+800)>size(dat,2)
        tmp(1:1001,i)=nan;
      elseif (saccs(i,1)-200)<1
        tmp(1:1001,i)=nan;
      else
        k = k +1;
        tmp(:,i) = pupil(saccs(i,1)-200:saccs(i,1)+800);
        meg(:,:,k) = dat(:,saccs(i,1)-200:saccs(i,1)+800);
      end
    end
    pupil_locked=nanmean(tmp,2);
    
    % COMPUTE FFT
    % ------------
    segleng = 200;
    segshift = 20;
    nseg=floor((size(meg,2)-segleng)/segshift+1);
    t = -100:20:1001;
    pxx = zeros(101,size(meg,1),nseg,'single');
    tmp_csd = zeros(size(meg,1),size(meg,1),101);
    pupil_count = 0; csd_count = 0; pxx_counter = 0;
    pup_seg = zeros(41,1);
    for isacc = 1 : size(meg,3)
      pupil_count = pupil_count+1; 
      isacc
      tmp_pxx = nan(101,size(meg,1),nseg);
      for iseg = 1 : nseg
        meg_seg  = squeeze(meg(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng,isacc));
        win = (hanning(200)./sum(hanning(200)));
        tmp_pup = squeeze(pupil((iseg-1)*segshift+1:(iseg-1)*segshift+segleng));
        pup_seg(iseg)  = pup_seg(iseg)+sum(win.*tmp_pup);
        
        if any(isnan(squeeze(meg_seg(1,:))))
          continue
        else
          [tmp_pxx(:,:,iseg),fxx] = pwelch(meg_seg',hanning(200),[],200,400,'power'); 
        end
      end
      if any(isnan(tmp_pxx(1,1,:)))
        continue
      else
        pxx_counter=pxx_counter+1;
        pxx = pxx + tmp_pxx;
      end
    end
    
    pup_seg = pup_seg ./ pupil_count; 
    
    for ifreq = 1 : size(pxx,1)
    	csd(:,:,ifreq) = sqrt(squeeze(pxx(ifreq,:,:)))*sqrt(squeeze(pxx(ifreq,:,:)))';
    end
  
%     pxx=nanmean(pxx,4);
%     pxx=(pxx-pxx(:,:,1))./pxx(:,:,1);
%     all_pxx(:,:,subj_counter) = squeeze(nanmean(pxx,2));
    

    n_saccs(subj_counter) = size(saccs,1);
    
    save([outdir sprintf('pp_revision_saccades_TFR_gla_isubj%d_iblock%d.mat',isubj,iblock)],'pup_seg','fxx','pupil_locked','pxx','-v7.3')
    clear pxx_src csd nai_src  pup_seg csd

  end
end

error('!')
%%
subj_counter = 10;
for isubj = 1:24
  if any(isubj==[2,3,5,6,21,9])
%   if any(isubj==[5,9])
    continue
  else
    subj_counter = subj_counter + 1;
  end
  load([outdir sprintf('pp_revision_saccades_TFR_gla_isubj%d_iblock%d.mat',isubj,iblock)])

  pxx_all(:,:,subj_counter,iblock)=100*(squeeze(nanmean(pxx(2:end,:,:),2))-squeeze(nanmean(pxx(2:end,:,1),2)))./squeeze(nanmean(pxx(2:end,:,1),2));

end


figure_w;

imagesc(nanmean(nanmean(pxx_all,3),4),[-2.5 2.5]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 101],'color','k','linestyle','--')
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:100,'yticklabel',fxx(2:5:end))
colormap(redblue); axis([1 41 1 101])
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_TFR_hh.pdf'))


%%

figure_w; hold on
% all_tmp = zscore(all_tmp);
all_tmp = all_tmp-nanmean(all_tmp(1:100,:),1);
s = std(all_tmp,[],2)/sqrt(size(all_tmp,2));
shadedErrorBar(1:901,nanmean(all_tmp,2),s)
line([0 901],[0 0],'color','k','linestyle',':')
line([100 100],[-2 2],'color','k','linestyle',':')
xlabel('Time [ms]'); ylabel('Pupil size')
tp_editplots; axis([0 800 -2 2])
set(gcf,'color','w')
set(gca,'xtick',[100 300 500 700 901],'xticklabel',[0 500 1000 1500 2000])

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_glasgow.pdf'))

conc_dat = cat(2,conc_dat,all_tmp);


%%


clear all_tmp pup_seg all_pxx

d=dir('~/pp/data_gla/fw4bt/osfstorage/data/ms01/pupil/*pupil_preproc_lp4.mat');

SUBJLIST = [];
for i = 1 : length(d)
  SUBJLIST = [SUBJLIST; d(i).name(1:4)];
end 


for isubj = 1:size(SUBJLIST,1)
  
  clear data dat pupil pup dataf src_r 
  
  for iblock = 1:1

    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load pupil data
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/pupil/%s_pupil_preproc_lp4.mat',SUBJLIST(isubj,:)))
      pupil(:,4) = data.trial{1}(4,:)';
      pupil(:,3) = data.trial{1}(3,:)';
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/pupil/%s_pupilgaze_preproc_lp4.mat',SUBJLIST(isubj,:)))
      pupil(:,1:2) = data.trial{1}(1:2,:)';
      %       pupil(:,1:2) =
    catch me
      continue
    end
    
    load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/cleanmeg_%s.mat',SUBJLIST(isubj,:)))
    data = cleanmeg; clear cleanmeg
    cfg = [];
    cfg.resamplefs = 400;
    data = ft_resampledata(cfg,data);
    f_sample = data.fsample;
            
    dat = data.trial{1};  clear data

    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    tmp = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil(:,4)=tmp;
    pupil = resample(pupil,400,600);
    
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------
    
    saccs=tp_detect_microsaccades(pupil(:,1:2),f_sample,5,6);
    
    rH(isubj) = corr(pupil(:,1),pupil(:,4));
    rV(isubj) = corr(pupil(:,2),pupil(:,4));
    
    [n,k]=hist(pupil(:,1),200); [~,i]=max(n); fix(1,isubj) = k(i);
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(2,isubj) = k(i);
    
    %     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
    dist = sqrt(((pupil(:,1)-fix(1,isubj)).^2) + ((pupil(:,2)-fix(2,isubj)).^2) );
    
    rExc(isubj) = corr(dist,pupil(:,4));
    pupil = zscore(pupil(:,4));
    
    tmp = []; meg = []; k = 0;
    for i = 1 : size(saccs,1)
      if (saccs(i,1)+800)>size(pupil,1)
        tmp(1:1001,i)=nan;
        continue
      elseif (saccs(i,1)-200)<1
        tmp(1:1001,i)=nan;
        continue
      elseif (saccs(i,1)+800)>size(dat,2)
        tmp(1:1001,i)=nan;
        continue
      else
        k = k +1;
        tmp(:,i) = pupil(saccs(i,1)-200:saccs(i,1)+800);
        meg(:,:,k) = dat(:,saccs(i,1)-200:saccs(i,1)+800);
      end
    end
    all_tmp(:,isubj,iblock)=nanmean(tmp,2);
    
    % COMPUTE FFT
    % ------------
    segleng = 200;
    segshift = 20;
    nseg=floor((size(meg,2)-segleng)/segshift+1);
    t = -100:100:1900
    pxx = nan(101,size(meg,1),nseg,size(meg,3),'single');
    for isacc = 1 : size(meg,3)
      isacc
      for iseg = 1 : nseg
        meg_seg  = squeeze(meg(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng,isacc));
        pup_seg(:,iseg,isubj,iblock)  = squeeze(pupil((iseg-1)*segshift+1:(iseg-1)*segshift+segleng));

        if any(isnan(squeeze(meg_seg(1,:))))
          pxx(:,:,iseg,isacc) = nan(101,size(meg_seg,1));
          continue
        else
          [pxx(:,:,iseg,isacc),fxx] = pwelch(meg_seg',hanning(200),[],200,400,'power');
        end
      end
    end
    pxx=nanmean(pxx,4);
    pxx=(pxx-pxx(:,:,1))./pxx(:,:,1);
    all_pxx(:,:,isubj) = squeeze(nanmean(pxx,2));
    
    n_saccs(isubj) = size(saccs,1);
    
  end
end
% 
% figure_w; hold on
% % all_tmp = zscore(all_tmp);
% all_tmp = all_tmp-nanmean(all_tmp(1:100,:),1);
% s = std(all_tmp,[],2)/sqrt(size(all_tmp,2));
% shadedErrorBar(1:901,nanmean(all_tmp,2),s)
% line([0 901],[0 0],'color','k','linestyle',':')
% line([100 100],[-2 2],'color','k','linestyle',':')
% xlabel('Time [ms]'); ylabel('Pupil size')
% tp_editplots; axis([0 800 -2 2])
% set(gcf,'color','w')
% set(gca,'xtick',[100 300 500 700 901],'xticklabel',[0 500 1000 1500 2000])
% 
% print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_muenster.pdf'))
% 
% conc_dat = cat(2,conc_dat,all_tmp);
% 
% figure_w; hold on
% % all_tmp = 100*(all_tmp-nanmean(all_tmp(1:100,:),1))./nanmean(all_tmp(1:100,:),1);
% s = std(conc_dat,[],2)/sqrt(size(conc_dat,2));
% shadedErrorBar(1:901,nanmean(conc_dat,2),s)
% line([0 901],[0 0],'color','k','linestyle',':')
% line([100 100],[-2 2],'color','k','linestyle',':')
% xlabel('Time [ms]'); ylabel('Pupil size')
% tp_editplots; axis([0 800 -.1 .1])
% set(gcf,'color','w')
% set(gca,'xtick',[100 300 500 700 901],'xticklabel',[0 500 1000 1500 2000])
% 
% print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_all.pdf'))

save([outdir 'pp_revision_saccades_TFR_mue.mat'],'pup_seg','fxx','all_pxx','all_tmp')

%%
load([outdir 'pp_revision_saccades_TFR_mue.mat'])

figure_w;

imagesc(nanmean(nanmean(all_pxx(:,:,:,:),3),4),[-0.025 0.025]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 101],'color','k','linestyle','--')
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:101,'yticklabel',fxx(1:5:end))
colormap(redblue)
% conc_dat = all_tmp;

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_TFR_hh.pdf'))

pup_seg = all_tmp;

pxx = nanmean(nanmean(all_pxx(:,:,:,:),3),4);