%% pp_mue_src_powerspectra

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
% v = 1;
% lag = 0;
% win_len = 800;
% overlap = 2; % 50% overlap
% -------------------------
% VERSION 2: with pupil lag
% -------------------------
% v = 2;
% lag = 1;
% win_len = 800;
% overlap = 2; % 50% overlap
% -------------------------
% VERSION 11: no pupil lag, less overlap
% -------------------------
% v = 11;
% lag = 0;
% win_len = 800;
% overlap = 1; % 0% overlap
% -------------------------
% VERSION 2: with pupil lag
% -------------------------
v = 22;
lag = 1;
win_len = 800;
overlap = 1; % 0% overlap
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

addpath /home/gnolte/meth/highlevel/

d=dir('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/*mat');

SUBJLIST = [];
for i = 1 : length(d)
  SUBJLIST = [SUBJLIST; d(i).name(end-7:end-4)];
end

addpath ~/pp/matlab/
trans = pp_transfer_gla2hh;
freqoi=2.^(1:(1/4):7);

%%
% -------------------------
for isubj =1:size(SUBJLIST,1)
  fprintf('Processing subj%d ...\n',isubj);
  
  clear pxx fxx pup pup_df
  
  for iblock = 1:1
    %
    fn = sprintf('pp_mue_src_powerspectra_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load pupil data
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/pupil/rawpupil_%s.mat',SUBJLIST(isubj,:)))
      
      pupil = puptc(:);
      
      % load meg data
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/cleanmeg_%s.mat',SUBJLIST(isubj,:)))
      data = cleanmeg; clear cleanmeg
      
      cfg=[];
      cfg.layout='CTF275.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:275),data.label);
      
    catch me
      %       src_r = nan(246,25);
      %       save([outdir fn '.mat'],'src_r')
      continue
    end

    k = 2;
    fnq = data.fsample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil);
    pupil = resample(pupil,400,600);
    
    
    data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    data.avg = data.trial{1}'; data.trial{1} = [];
    
    data.avg = resample(data.avg,400,600);
    f_sample = 400;
    
    if lag
      pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
      pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    
    pupil_df = diff(pupil);
    
    
    load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/leadfields/lf_%s.mat',SUBJLIST(isubj,:)))
    
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)
      
      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = f_sample;
      para.overlap = 0.5;
      csd(:,:,ifreq)=tp_compute_csd_wavelets(data.avg',para);
      
    end
    
    csd = nanmean(csd,3);
    
    % -------------------------------
    % beamforming
    % -------------------------------
    para      = [];
    para.iscs = 1;
    para.reg  = 0.05;
    filt   = tp_beamformer(real(csd),lf,para);
    % -------------------------------
    
    opt.n_win = win_len; % 10s segment length, i.e., 0.1:0.1:100
    opt.n_shift = win_len/overlap; % no overlap
    
    nseg=floor((size(data.avg,1)-opt.n_win)/opt.n_shift+1);
    clear pxx fxx pup pup_df
    ff = 2:1/(opt.n_win/400):128;
    
    pupil = pupil(1:size(data.avg,1));
    pup_nanidx = isnan(pupil);
    pupil_df = diff(pupil);
    
    data.avg(pup_nanidx,:)=nan;
    
    pxx = nan(size(ff,2),max(BNA.tissue_5mm(:)),nseg);
    for iseg = 1 : nseg
      fprintf('%d / %d\n',iseg,nseg)
      seg_dat = data.avg((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win,:)*filt;
      seg_dat = seg_dat(:,trans);
      
      if any(isnan(seg_dat(:,1)))
        pup(iseg) = nan;
        pup_df(iseg)=nan;
        continue
      end
      
      [tmp_pxx,fxx]=pwelch(seg_dat,hanning(opt.n_win),0.5,ff,400,'power');
      
      for igrid = 1 : max(BNA.tissue_5mm(:))
        pxx(:,igrid,iseg) = mean(tmp_pxx(:,BNA.tissue_5mm == igrid),2);
      end
      
      pup(iseg)  = mean(pupil((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
      if iseg ~= nseg
        pup_df(iseg) = mean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
      else
        pup_df(iseg) = mean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win-1));
      end
    end
    
    save([outdir fn '.mat'],'pxx','fxx','pup','pup_df')
    
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

error('!')
