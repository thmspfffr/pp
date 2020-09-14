%% pp_gla_src_scaling
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
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_gla_src_scaling_s%d_b%d_v%d',isubj,iblock,v);
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
%       
%       [outp.pxx,outp.fxx]=pwelch(data.trial{1}',hanning(400),0,1:1:200,400);
  
    catch me
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
    
    pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
    pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    
%     data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    
%     tmp = data;
    data.avg = data.trial{1}'; %data.trial{1} = [];

    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
   
                
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)
    
      for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
      end

      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;  
      para.meth     = 'conv';
      [tp_csd, ~,opt,dataf]=tp_compute_csd_wavelets(data.avg',para);
      
      % -------------------------------
      % compute power from csd
      % -------------------------------
      tmp = diag(abs(tp_csd)); 
      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      pup = pupil;
      % -------------------------------
      % beamforming
      % -------------------------------
      para      = [];
      para.iscs = 1;
      para.reg  = 0.05;
      [tp_filt,tp_pow]      = tp_beamformer(real(tp_csd),lf,para);
      % -------------------------------
%             
      n_win = 20 * 400;
      n_shift = n_win/2;
      nseg=floor((size(dataf,2)-n_win)/n_shift+1);

      for iwin = 1 : nseg
                
        iwin
        
        tmp1= dataf(:,(iwin-1)*n_shift+1:(iwin-1)*n_shift+n_win);
        tmp1 = abs(tp_filt'*tmp1).^2;
        tmp2 = pup((iwin-1)*n_shift+1:(iwin-1)*n_shift+n_win);
        
        if any(isnan(tmp1(1,:))) || any(isnan(tmp2))
          idx(iwin)=0;
          continue
        else 
          idx(iwin)=1;
        end
          
        
        [pxx,fxx]=pwelch(tmp1',hanning(n_win),0,0.05:0.05:2,400);
        X = [ones(1,length(fxx))' log10(fxx)];
        Y = log10(pxx);
        tmp = X\Y; slope_meg(:,iwin) = tmp(2,:)';
        
        [pxx,fxx]=pwelch(tmp2',hanning(n_win),0,0.05:0.05:2,400);
        X = [ones(1,length(fxx))' log10(fxx')];
        Y = log10(pxx);
        tmp = X\Y'; slope_pup(:,iwin) = tmp(2);
        
  
      end
      
      outp.corr(ifreq)=corr(slope_meg(:,idx)',slope_pup(idx)');  
      clear src pup csd tp_csd dataf 
      
    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

error('!')
%%
