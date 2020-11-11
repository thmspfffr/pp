%% pp_gla_src_fooof
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
v = 1;
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
win_len = 1600;
lag = 0;
% -------------------------
% VERSION 4: with pupil lag
% -------------------------
% v = 2;
% SUBJLIST  = 1:24;
% freqoi    = 2.^(1:(1/4):7);
% win_len = 1600;
% lag = 1;
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
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_gla_src_fooof_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
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
      
      chanidx = outp.chanidx;
      save(sprintf('~/pp/proc/src/chanidx_s%d.mat',isubj),'chanidx')
      
%       
      [outp.pxx,outp.fxx]=pwelch(data.trial{1}',hanning(400),0,1:1:200,400);
            
      outp.tp_sens_pow = nan(248,25);
      outp.sens_r      = nan(248,25);
      outp.sens_r_sp   = nan(248,25);
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
%     len = min([size(data.trial{1},2) size(pupil,1)]);
%     data.trial{1} = data.trial{1}(:,1:len);
%     pupil = pupil(1:len,:);
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    if lag
        pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
        pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    
%     data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    
%     tmp = data;
    data.avg = data.trial{1}'; %data.trial{1} = [];
%     clear data

    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
    
    
    for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
    end

    clear tp_csd
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)

      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;  
      para.overlap = 0.5;
      tp_csd(:,:,ifreq)=tp_compute_csd_wavelets(data.avg',para);
      
    end
    
    tp_csd = nanmean(tp_csd,3);
      
    % -------------------------------
    % beamforming
    % -------------------------------
    para      = [];
    para.iscs = 1;
    para.reg  = 0.05;
    tp_filt   = tp_beamformer(real(tp_csd),lf,para);
    % -------------------------------
    
    data_src = data.avg*tp_filt; 
    clear data
    
    opt.n_win = win_len; % 10s segment length, i.e., 0.1:0.1:100
    opt.n_shift = win_len; % no overlap
    
    nseg=floor((size(data_src,1)-opt.n_win)/opt.n_shift+1);
    clear pxx fxx pup pup_df
    ff = 3:1/(opt.n_win/400):50;
    
    pupil = pupil(1:size(data_src,1));
    pup_nanidx = isnan(pupil);
    pupil_df = diff(pupil);
    
    data_src(pup_nanidx,:)=nan;
     
    pxx = nan(size(ff,2),size(tp_filt,2),nseg);
    for iseg = 1 : nseg
        fprintf('%d / %d\n',iseg,nseg)
        seg_dat = data_src((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win,:);
        
        if any(isnan(seg_dat(:,1)))
%             pxx(:,:,iseg) = nan(size(fxx,1),size(tp_filt,2));
            pup(iseg) = nan;
            pup_df(iseg)=nan;
            continue        
        end
        
        [pxx(:,:,iseg),fxx]=pwelch(seg_dat,hanning(opt.n_win),[],ff,400,'power');
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
%%
