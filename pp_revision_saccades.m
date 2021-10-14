%% pp_src_pupil_power_correlations
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
for isubj = 25:34
  
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    %
    fn = sprintf('pp_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
%     if tp_parallel(fn,outdir,1,0)
%       continue
%     end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
      % data from 'pp_prepare_data.m'
        load(sprintf('~/pp/data/ham/pp_rest_s%d_b%d_v%d.mat',isubj,iblock,1))
%       load(sprintf('~/pp/data/ham/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
    catch me
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
    
    tmp = filtfilt(bhil, ahil, pupil(:,[1 end]));
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
    % ------
    saccs=tp_detect_microsaccades(pupil(:,2:3),400,5);
    
    rH(isubj) = corr(pupil(:,2),pupil(:,4));
    rV(isubj) = corr(pupil(:,3),pupil(:,4));
    
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(1,isubj) = k(i);
    [n,k]=hist(pupil(:,3),200); [~,i]=max(n); fix(2,isubj) = k(i);

%     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
    dist = sqrt(((pupil(:,2)-fix(1,isubj)).^2) + ((pupil(:,3)-fix(2,isubj)).^2) );
    
    rExc(isubj,iblock) = corr(dist(~isnan(pupil(:,3))),pupil(~isnan(pupil(:,3)),4));

    
    n_saccs(isubj,iblock) = size(saccs,1);
    
        pupil = pupil(:,end);

%     loc_pup = [];  sacc_rate = [];
%     segleng = 20*400;
%     segshift = segleng;
%     
%     nseg=floor((size(pupil,1)-segleng)/segshift+1);
%         
%     for iseg = 1 : nseg
% 
%       loc_pup(iseg)=mean(pupil((iseg-1)*segleng+1:(iseg-1)*segshift+segleng));
%       sacc_rate(iseg) = sum(saccs(:,1)>=(iseg-1)*segleng+1 & saccs(:,1)<= (iseg-1)*segshift+segleng);
% 
%     end
    tmp = []
    for i = 7 : size(saccs,1)
      if (saccs(i,1)+800)>240000
        continue
      end
      tmp(:,i-6) = pupil(saccs(i,1)-100:saccs(i,1)+800);
    end
    
    all_tmp(:,isubj,iblock)=mean(tmp,2);
 
  end
  
end

figure; hold on
s = std(nanmean(all_tmp(:,25:34,:),3),[],2);
shadedErrorBar(1:901,mean(nanmean(all_tmp(:,25:34,:),3),2),s)
line([0 901],[0 0],'color','k','linestyle','--')


%% GLASGOW DATA
% -------------------------



for isubj = 1:24
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
%     fn = sprintf('pp_gla_src_powerspectra_s%d_b%d_v%d',isubj,iblock,v);
%     if tp_parallel(fn,outdir,1,0)
%       continue
%     end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
%     try
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pupil = data.trial{1}';
      f_sample = data.fsample;
     
%     catch me
%       src_r = nan(246,25);
%       save([outdir fn '.mat'],'src_r')
%       continue
%     end
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    tmp = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil(:,4)=tmp;
    pupil = resample(pupil,400,1000);
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------
    
    saccs=tp_detect_microsaccades(pupil(:,1:2),f_sample,5);

    rH(isubj) = corr(pupil(:,1),pupil(:,4));
    rV(isubj) = corr(pupil(:,2),pupil(:,4));
    
    [n,k]=hist(pupil(:,1),200); [~,i]=max(n); fix(1,isubj) = k(i);
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(2,isubj) = k(i);

%     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
    dist = sqrt(((pupil(:,1)-fix(1,isubj)).^2) + ((pupil(:,2)-fix(2,isubj)).^2) );
    
    rExc(isubj) = corr(dist,pupil(:,4));
    
    tmp = []
    for i = 4 : size(saccs,1)
      if (saccs(i,1)+800)>180000
        continue
      end
      tmp(:,i-3) = pupil(saccs(i,1)-100:saccs(i,1)+800);
    end
    
    all_tmp(:,isubj)=mean(tmp,2);
    
    n_saccs(isubj) = size(saccs,1);

  end
end

% if isubj ~= 34
%     continue
%     else
      figure; hold on
      s = std(nanmean(all_tmp(:,1:end,:),3),[],2);
      shadedErrorBar(1:901,mean(nanmean(all_tmp(:,1:end,:),3),2),s)
      line([0 1601],[0 0],'color','k','linestyle','--')
%       
%       
%     end


%%

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

d=dir('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/*mat');

SUBJLIST = [];
for i = 1 : length(d)
  SUBJLIST = [SUBJLIST; d(i).name(end-7:end-4)];
end


for isubj = 1:size(SUBJLIST,1)
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
%     fn = sprintf('pp_gla_src_powerspectra_s%d_b%d_v%d',isubj,iblock,v);
%     if tp_parallel(fn,outdir,1,0)
%       continue
%     end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
%     try
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pupil = data.trial{1}';
      f_sample = data.fsample;
     
%     catch me
%       src_r = nan(246,25);
%       save([outdir fn '.mat'],'src_r')
%       continue
%     end
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    tmp = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil(:,4)=tmp;
    pupil = resample(pupil,400,1000);
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------
    
    saccs=tp_detect_microsaccades(pupil(:,1:2),f_sample,5);

    rH(isubj) = corr(pupil(:,1),pupil(:,4));
    rV(isubj) = corr(pupil(:,2),pupil(:,4));
    
    [n,k]=hist(pupil(:,1),200); [~,i]=max(n); fix(1,isubj) = k(i);
    [n,k]=hist(pupil(:,2),200); [~,i]=max(n); fix(2,isubj) = k(i);

%     fix(:,isubj) = [nanmean(pupil(:,1)),nanmean(pupil(:,2))];
    dist = sqrt(((pupil(:,1)-fix(1,isubj)).^2) + ((pupil(:,2)-fix(2,isubj)).^2) );
    
    rExc(isubj) = corr(dist,pupil(:,4));
    
    tmp = []
    for i = 4 : size(saccs,1)
      tmp(:,i-1) = pupil(saccs(i,1)-100:saccs(i,1)+800);
    end
    
    all_tmp(:,isubj)=mean(tmp,2);

  end
end

% if isubj ~= 34
%     continue
%     else
      figure; hold on
      s = std(nanmean(all_tmp(:,1:end,:),3),[],2);
      shadedErrorBar(1:901,mean(nanmean(all_tmp(:,1:end,:),3),2),s)
      line([0 1601],[0 0],'color','k','linestyle','--')
%       
%       
%     end

