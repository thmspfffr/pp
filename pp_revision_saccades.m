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
clear rExc
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
    catch me
      disp('Sth went wrong')
      percent_rejected(subj_counter,iblock) =nan;
      continue
    end
        
    cfg=[];
    cfg.layout='CTF275.lay';
    lay = ft_prepare_layout(cfg);
    [~, outp.chanidx] = ismember(lay.label(1:275),label(startsWith(label,'M')));
%     size(dat)

    idx=zeros(size(pupil,1),1);
    for iart = 1 : size(art.blinks_ts,1)
      if (art.blinks_ts(iart,1)-200)<=0
        idx(1:art.blinks_ts(iart,1)+200)=1;
      elseif art.blinks_ts(iart,2)+200>size(pupil,1)
        idx(art.blinks_ts(iart,1):end)=1;
      else
        idx(art.blinks_ts(iart,1)-200:art.blinks_ts(iart,2)+200)=1;
      end
    end
    pupil(:,6)=idx;

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
   
    percent_rejected(subj_counter,iblock) = 100*sum(pupil(:,6))/size(pupil,1);
%     % ------
    saccs=tp_detect_microsaccades(pupil(:,2:3),400,5,6);
        
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(1,subj_counter) = k(i);
    [n,k]=hist(pupil(:,3),200); [~,i]=max(n); fix(2,subj_counter) = k(i);
%     
    dist = sqrt(((pupil(:,2)-fix(1,subj_counter)).^2) + ((pupil(:,3)-fix(2,subj_counter)).^2) );
%     
    rExc(subj_counter,iblock) = corr(dist(~isnan(pupil(:,3))),pupil(~isnan(pupil(:,3)),4));
%     
    pupil1 = zscore(pupil(:,5));
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
        tmp1(:,i) = pupil1(saccs(i,1)-200:saccs(i,1)+800);
        meg(:,:,k) = dat(:,saccs(i,1)-200:saccs(i,1)+800);
      end
    end
    pupil_locked=nanmean(tmp,2);
    pupil_locked_noregression=nanmean(tmp1,2);
    
%     COMPUTE FFT
%     ------------
    segleng = 200;
    segshift = 20;
    nseg=floor((size(meg,2)-segleng)/segshift+1);
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

    save([outdir sprintf('pp_revision_saccades_TFR_hh_isubj%d_iblock%d.mat',isubj,iblock)],'pup_seg','fxx','pupil_locked','pupil_locked_noregression','pxx','-v7.3')

    clear pxx_src csd nai_src  pup_seg csd
    
  end
end
rExc = mean(rExc,2);
save('~/pp/proc/src/pp_hh_rexc.mat','rExc')

error('!')

%% GLASGOW DATA
% -------------------------
% subj2,3,5, 21: impossible
    % subj6,, doable?
rExc = []; clear all_tmp all_pxx n_saccs percent_rejected rExc fix

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
    
    
%     saccs = data.cfg.artfctdef.sacc.artifact;
    pupil = data.trial{1}';
    f_sample = data.fsample;
    
    art.blinks_ts=data.cfg.artfctdef.blink.artifact;
    idx=zeros(size(pupil,1),1);
    for iart = 1 : size(art.blinks_ts,1)
      if (art.blinks_ts(iart,1)-200)<=0
        idx(1:art.blinks_ts(iart,1)+200)=1;
      elseif art.blinks_ts(iart,2)+200>size(pupil,1)
        idx(art.blinks_ts(iart,1):end)=1;
      else
        idx(art.blinks_ts(iart,1)-200:art.blinks_ts(iart,2)+200)=1;
      end
    end
    pupil(:,6)=idx;
    
    size(art.blinks_ts)
    
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
    
    [n,k]=hist(pupil(:,1),200); [~,i]=max(n); fix(1,subj_counter) = k(i);
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(2,subj_counter) = k(i);
    figure_w; hold on
    plot(pupil(:,1),pupil(:,2),'.');
    plot(fix(1,subj_counter),fix(2,subj_counter),'r.','markersize',30)
    dist = sqrt(((pupil(:,1)-fix(1,subj_counter)).^2) + ((pupil(:,2)-fix(2,subj_counter)).^2) );
    
    rExc(subj_counter) = corr(dist,pupil(:,4));
    percent_rejected(subj_counter) = 100*sum(pupil(:,6))/size(pupil,1);

    pupil = zscore(pupil(:,4));
%     
%     tmp = []; meg = []; k = 0;
%     for i = 1 : size(saccs,1)
%       if (saccs(i,1)+800)>size(dat,2)
%         tmp(1:1001,i)=nan;
%       elseif (saccs(i,1)-200)<1
%         tmp(1:1001,i)=nan;
%       else
%         k = k +1;
%         tmp(:,i) = pupil(saccs(i,1)-200:saccs(i,1)+800);
%         meg(:,:,k) = dat(:,saccs(i,1)-200:saccs(i,1)+800);
%       end
%     end
%     pupil_locked=nanmean(tmp,2);
%     
%     % COMPUTE FFT
%     % ------------
%     segleng = 200;
%     segshift = 20;
%     nseg=floor((size(meg,2)-segleng)/segshift+1);
%     t = -100:20:1001;
%     pxx = zeros(101,size(meg,1),nseg,'single');
%     tmp_csd = zeros(size(meg,1),size(meg,1),101);
%     pupil_count = 0; csd_count = 0; pxx_counter = 0;
%     pup_seg = zeros(41,1);
%     for isacc = 1 : size(meg,3)
%       pupil_count = pupil_count+1; 
%       isacc
%       tmp_pxx = nan(101,size(meg,1),nseg);
%       for iseg = 1 : nseg
%         meg_seg  = squeeze(meg(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng,isacc));
%         win = (hanning(200)./sum(hanning(200)));
%         tmp_pup = squeeze(pupil((iseg-1)*segshift+1:(iseg-1)*segshift+segleng));
%         pup_seg(iseg)  = pup_seg(iseg)+sum(win.*tmp_pup);
%         
%         if any(isnan(squeeze(meg_seg(1,:))))
%           continue
%         else
%           [tmp_pxx(:,:,iseg),fxx] = pwelch(meg_seg',hanning(200),[],200,400,'power'); 
%         end
%       end
%       if any(isnan(tmp_pxx(1,1,:)))
%         continue
%       else
%         pxx_counter=pxx_counter+1;
%         pxx = pxx + tmp_pxx;
%       end
%     end
%     
%     pup_seg = pup_seg ./ pupil_count; 
%     
%     for ifreq = 1 : size(pxx,1)
%     	csd(:,:,ifreq) = sqrt(squeeze(pxx(ifreq,:,:)))*sqrt(squeeze(pxx(ifreq,:,:)))';
%     end
%  
%     n_saccs(subj_counter) = size(saccs,1);
%     
%     save([outdir sprintf('pp_revision_saccades_TFR_gla_isubj%d_iblock%d.mat',isubj,iblock)],'pup_seg','fxx','pupil_locked','pxx','-v7.3')
    clear pxx_src csd nai_src  pup_seg csd

  end
end

save('~/pp/proc/src/pp_gla_rexc.mat','rExc')

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

d=dir('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/*mat');

SUBJLIST = [];
for i = 1 : length(d) 
    SUBJLIST = [SUBJLIST; d(i).name(end-7:end-4)];  
end

rExc = []; clear all_tmp all_pxx n_saccs percent_rejected fix rExc

subj_counter = 0;  

for isubj = 1:size(SUBJLIST,1)
  subj_counter = subj_counter + 1;
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
    
    
    art.blinks_ts=data.cfg.artfctdef.blink.artifact;
    idx=zeros(size(pupil,1),1);
    for iart = 1 : size(art.blinks_ts,1)
      if (art.blinks_ts(iart,1)-200)<=0
        idx(1:art.blinks_ts(iart,1)+200)=1;
      elseif art.blinks_ts(iart,2)+200>size(pupil,1)
        idx(art.blinks_ts(iart,1):end)=1;
      else
        idx(art.blinks_ts(iart,1)-200:art.blinks_ts(iart,2)+200)=1;
      end
    end
    pupil(:,6)=idx;
    
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
       
    [n,k]=hist(pupil(:,1),200); [~,i]=max(n); fix(1,subj_counter) = k(i);
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(2,subj_counter) = k(i);
    
    figure_w; hold on
    plot(pupil(:,1),pupil(:,2),'.');
    plot(fix(1,subj_counter),fix(2,subj_counter),'r.','markersize',30)
   
%     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
    dist = sqrt(((pupil(:,1)-fix(1,subj_counter)).^2) + ((pupil(:,2)-fix(2,subj_counter)).^2) );
    
    rExc(subj_counter) = corr(dist,pupil(:,4));
    
        percent_rejected(subj_counter,iblock) = 100*sum(pupil(:,6))/size(pupil,1);
% 
    pupil = zscore(pupil(:,4));
    
%     tmp = []; meg = []; k = 0;
%     for i = 1 : size(saccs,1)
%       if (saccs(i,1)+800)>size(dat,2)
%         tmp(1:1001,i)=nan;
%       elseif (saccs(i,1)-200)<1
%         tmp(1:1001,i)=nan;
%       else
%         k = k +1;
%         tmp(:,i) = pupil(saccs(i,1)-200:saccs(i,1)+800);
%         meg(:,:,k) = dat(:,saccs(i,1)-200:saccs(i,1)+800);
%       end
%     end
%     pupil_locked=nanmean(tmp,2);
    
%     % COMPUTE FFT
%     % ------------
%     segleng = 200;
%     segshift = 20;
%     nseg=floor((size(meg,2)-segleng)/segshift+1);
%     t = -100:20:1001;
%     pxx = zeros(101,size(meg,1),nseg,'single');
%     tmp_csd = zeros(size(meg,1),size(meg,1),101);
%     pupil_count = 0; csd_count = 0; pxx_counter = 0;
%     pup_seg = zeros(41,1);
%     for isacc = 1 : size(meg,3)
%       pupil_count = pupil_count+1; 
%       isacc
%       tmp_pxx = nan(101,size(meg,1),nseg);
%       for iseg = 1 : nseg
%         meg_seg  = squeeze(meg(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng,isacc));
%         win = (hanning(200)./sum(hanning(200)));
%         tmp_pup = squeeze(pupil((iseg-1)*segshift+1:(iseg-1)*segshift+segleng));
%         pup_seg(iseg)  = pup_seg(iseg)+sum(win.*tmp_pup);
%         
%         if any(isnan(squeeze(meg_seg(1,:))))
%           continue
%         else
%           [tmp_pxx(:,:,iseg),fxx] = pwelch(meg_seg',hanning(200),[],200,400,'power'); 
%         end
%       end
%       if any(isnan(tmp_pxx(1,1,:)))
%         continue
%       else
%         pxx_counter=pxx_counter+1;
%         pxx = pxx + tmp_pxx;
%       end
%     end
%     
%     pup_seg = pup_seg ./ pupil_count; 
%     
%     for ifreq = 1 : size(pxx,1)
%     	csd(:,:,ifreq) = sqrt(squeeze(pxx(ifreq,:,:)))*sqrt(squeeze(pxx(ifreq,:,:)))';
%     end
%     
%     save([outdir sprintf('pp_revision_saccades_TFR_mue_isubj%d_iblock%d.mat',isubj,iblock)],'pup_seg','fxx','pupil_locked','pxx','-v7.3')
    clear pxx_src csd nai_src  pup_seg csd
    
  end
end

save('~/pp/proc/src/pp_mue_rexc.mat','rExc')

%%

load '~/pp/proc/src/pp_gla_rexc.mat'
rrr = rExc';
load '~/pp/proc/src/pp_hh_rexc.mat'
rrr = [rrr; rExc];
load '~/pp/proc/src/pp_mue_rexc.mat'
rrr = [rrr; rExc'];

mean(rrr)
[~,p]=ttest(rrr)

%% PLOT MUENSTER
clear pxx_all
d=dir('~/pp/data_gla/fw4bt/osfstorage/data/ms01/pupil/*pupil_preproc_lp4.mat');

SUBJLIST = [];
for i = 1 : length(d)
  SUBJLIST = [SUBJLIST; d(i).name(1:4)];
end 

subj_counter = 0;  

for isubj = 1:size(SUBJLIST,1)
  subj_counter = subj_counter + 1;
  load([outdir sprintf('pp_revision_saccades_TFR_mue_isubj%d_iblock%d.mat',isubj,iblock)])
  pxx_all(:,:,subj_counter,iblock)=100*(squeeze(nanmean(pxx(2:end,:,:),2))-squeeze(nanmean(pxx(2:end,:,1),2)))./squeeze(nanmean(pxx(2:end,:,1),2));
end
pxx_all(:,:,23)=[];

%%
figure_w;
imagesc(nanmean(nanmean(pxx_all(:,:,:),3),4),[-2.5 2.5]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 101],'color','k','linestyle','--')
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:100,'yticklabel',fxx(2:5:end))
colormap(redblue); axis([1 41 1 101])
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_TFR_hh.pdf'))


%% PLOT ALL DATA POOLED
% -----------------------------

clear rrr pxx_all pup_all

cfg=[];
cfg.layout='CTF275.lay';
lay = ft_prepare_layout(cfg);

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
subj_counter = 0; clear pxx_all pxx_src_all rrr pup_all
for isubj = 25:34
  subj_counter=subj_counter+1

  for iblock = 1 : 2
      load(sprintf('~/pp/data/ham/pp_rest_s%d_b%d_v%d.mat',isubj,iblock,1),'label')
      [~, outp.chanidx] = ismember(lay.label(1:275),label(startsWith(label,'M')));

      load([outdir sprintf('pp_revision_saccades_TFR_hh_isubj%d_iblock%d.mat',isubj,iblock)])
           
      tmp_pxx = nan(size(pxx,1)-1,275,size(pxx,3));
      tmp_pxx(:,outp.chanidx>0,:)=100*(pxx(2:end,:,:)-pxx(2:end,:,1))./pxx(2:end,:,1);
      
      minmax_hh=[min(lay.pos(1:275,2)) max(lay.pos(1:275,2))];
      ser_hh = linspace(minmax_hh(1),minmax_hh(2),40);

      for i = 1 : size(ser_hh,2)-1
        idx = lay.pos(1:275,2)>=ser_hh(i) & lay.pos(1:275,2)<ser_hh(i+1);
        pxx_all(:,:,i,subj_counter,iblock)      = squeeze(nanmean(tmp_pxx(:,idx,:),2));
      end
      
      for ifreq = 1 : 100
        rrr(:,ifreq,subj_counter,iblock) = corr(squeeze(pxx_all(ifreq,:,:,subj_counter,iblock)),pup_seg);
      end
      
      pup_all(:,subj_counter,iblock)=100*(pup_seg-pup_seg(1,:))./pup_seg(1,:);
  end
end

pxx_all = nanmean(pxx_all,5);
pup_all = nanmean(pup_all,3);
rrr = nanmean(rrr,4);

cfg=[];
cfg.layout='4D248.lay';
lay = ft_prepare_layout(cfg);
minmax_gla=[min(lay.pos(1:248,2)) max(lay.pos(1:248,2))];
ser_gla = linspace(minmax_gla(1),minmax_gla(2),40);

for isubj = 1:24
%   if any(isubj==[2,3,5,6,21,9])
  load(sprintf('~/pp/proc/src/chanidx_s%d.mat',isubj),'chanidx')

  if any(isubj==[5,9])
    continue
  else
    subj_counter = subj_counter + 1;
  end
  load([outdir sprintf('pp_revision_saccades_TFR_gla_isubj%d_iblock%d.mat',isubj,1)])

  tmp_pxx = nan(size(pxx,1)-1,248,size(pxx,3));
  tmp_pxx(:,chanidx>0,:)=100*(pxx(2:end,:,:)-pxx(2:end,:,1))./pxx(2:end,:,1);
      
  for i = 1 : size(ser_gla,2)-1
    idx = lay.pos(1:248,2)<ser_gla(i+1) & lay.pos(1:248,2)>ser_gla(i);
    pxx_all(:,:,i,subj_counter) = squeeze(nanmean(tmp_pxx(:,idx,:),2));
  end

  pup_all(:,subj_counter)=100*(pup_seg-pup_seg(1,:))./pup_seg(1,:);
  for ifreq = 1 : 100
    rrr(:,ifreq,subj_counter) = corr(squeeze(pxx_all(ifreq,:,:,subj_counter)),pup_seg);
  end
  
end

% LOAD MUENSTER
% -----------------------------------------------

% load 275 CTF layout
cfg=[];
cfg.layout='CTF275.lay';
lay = ft_prepare_layout(cfg);
SUBJLIST=1:41; SUBJLIST([10,12,17,19,22,27,35,38,39,40])=[];
% SUBJLIST=1:41; SUBJLIST(list)=[];

for n_subj = 1:length(SUBJLIST)
  isubj = SUBJLIST(n_subj);
  subj_counter = subj_counter + 1;
  load([outdir sprintf('pp_revision_saccades_TFR_mue_isubj%d_iblock%d.mat',isubj,1)])
  load(sprintf('~/pp/proc/src/chanidx_mue_s%d.mat',isubj),'chanidx')

  tmp_pxx = nan(size(pxx,1)-1,275,size(pxx,3));
  tmp_pxx(:,chanidx>0,:)=100*(pxx(2:end,:,:)-pxx(2:end,:,1))./pxx(2:end,:,1);
  
  minmax_hh=[min(lay.pos(1:275,2)) max(lay.pos(1:275,2))];
  ser_hh = linspace(minmax_hh(1),minmax_hh(2),40);
  
  for i = 1 : size(ser_hh,2)-1
    idx = lay.pos(1:275,2)>=ser_hh(i) & lay.pos(1:275,2)<ser_hh(i+1);
    pxx_all(:,:,i,subj_counter)      = squeeze(nanmean(tmp_pxx(:,idx,:),2));
  end
  
  for ifreq = 1 : 100
    rrr(:,ifreq,subj_counter) = corr(squeeze(pxx_all(ifreq,:,:,subj_counter)),pup_seg);
  end
  pup_all(:,subj_counter)=100*(pup_seg-pup_seg(1,:))./pup_seg(1,:);

end

pxx_all(:,:,:,end-12)=[];
pup_all(:,end-12) = [];
rrr(:,:,end-12)=[];
%%
t = -250:50:1750;

close 
figure_w; subplot(3,3,1);
imagesc(nanmean(nanmean(pxx_all(1:67,:,:,:),3),4),[-5 5]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 67],'color','k','linestyle',':')
line([13 13],[0 67],'color','k','linestyle','-')

set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:67,'yticklabel',fxx(2:5:end))
colormap(redblue); axis([1 41 1 67])

subplot(3,3,4)
imagesc(squeeze(nanmean(nanmean(pxx_all(1:67,13,:,:),2),4))',[-5 5]); set(gca,'ydir','normal')
set(gca,'xtick',1:5:67,'xticklabel',fxx(2:5:end))
colormap(redblue);
tp_editplots

[h,p]=ttest(squeeze(pxx_all(1:67,13,:,:)),zeros(size(squeeze(pxx_all(1:67,13,:,:)))),'dim',3); h=p<fdr1(p(:),0.1,0);

subplot(3,3,5)
imagesc(squeeze(nanmean(nanmean(pxx_all(1:67,13,:,:),2),4))'.*h',[-5 5]); set(gca,'ydir','normal')
set(gca,'xtick',1:5:67,'xticklabel',fxx(2:5:end))
colormap(redblue);
tp_editplots

% 
[h,p] = ttest(squeeze(nanmean(pxx_all,3)),zeros(size(squeeze(nanmean(pxx_all,3)))),'dim',3); h=p<fdr1(p(:),0.1,0);
subplot(3,3,2);
imagesc(nanmean(nanmean(pxx_all(1:67,:,:,:),3),4).*h(1:67,:,:),[-5 5]); set(gca,'ydir','normal')
tp_editplots
line([6 6],[0 62],'color','k','linestyle','--')
line([13 13],[0 67],'color','k','linestyle','-')
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
set(gca,'ytick',1:5:62,'yticklabel',fxx(2:5:end))
colormap(redblue); axis([1 41 1 62])


subplot(3,3,7);
imagesc(nanmean(rrr(:,1:67,:),3),[-0.5 0.5]); set(gca,'ydir','normal')
set(gca,'xtick',1:5:67,'xticklabel',fxx(2:5:end))
tp_editplots

subplot(3,3,8);
[h,p] = ttest(rrr, zeros(size(rrr)),'dim',3); h=p<fdr1(p(:),0.1,0);
imagesc(nanmean(rrr(:,1:67,:),3).*h(:,1:67),[-0.5 0.5]); set(gca,'ydir','normal')
set(gca,'xtick',1:5:67,'xticklabel',fxx(2:5:end))
tp_editplots


print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_TFR_all.pdf'))

g=figure_w;
[h,p]=ttest(pup_all,zeros(size(pup_all)),'dim',2);
errBar = nanstd(pup_all,[],2)./sqrt(size(pup_all,2));
shadedErrorBar([],nanmean(pup_all,2),errBar,'k',1)
line([6 6],[-300 300],'color','k','linestyle',':')
line([0 31],[0 0],'color','k','linestyle',':')
text([5],[70],sprintf('All P>%.3f',min(p)))
set(gca,'xtick',1:5:41,'xticklabel',t(1:5:end))
tp_editplots; axis([1 41 -300 300])
set(g,'Renderer','Painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_revision_saccades_pupil_all.pdf'))

