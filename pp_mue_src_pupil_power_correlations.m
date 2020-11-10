%% pp_mue_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
% v = 1;
% freqoi    = 2.^(1:(1/4):7);
% lag = 0;
% -------------------------
% VERSION 3: with pupil lag
% -------------------------
v = 2;
freqoi    = 2.^(1:(1/4):7);
lag = 1;
% -------------------------

addpath('~/Documents/MATLAB/fieldtrip-20181231/')
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

%%
% -------------------------
for isubj = 1:size(SUBJLIST,1)
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_mue_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
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
      f_sample = data.fsample;
   
      cfg=[];
      cfg.layout='CTF275.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:275),data.label);
    
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil);
    
    if lag
        pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
        pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    pupil_df = diff(pupil);
    
    data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    
%     tmp = data;
    data.avg = data.trial{1}'; %data.trial{1} = [];

    load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/leadfields/lf_%s.mat',SUBJLIST(isubj,:)))
    outp.sens_pow   = nan(275,25);
    outp.sens_r     = nan(275,25);
    outp.sens_r_df  = nan(275,25);
    
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)

      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = f_sample;  
      para.overlap = 0.5;
      [tp_csd, dataf,opt]=tp_compute_csd_wavelets(data.avg',para);
      
      % -------------------------------
      % compute power from csd
      % -------------------------------
      tmp = diag(abs(tp_csd)); 
      outp.sens_pow(outp.chanidx>0,ifreq) = tmp(outp.chanidx(outp.chanidx>0));
      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      nseg=floor((size(data.avg,1)-opt.n_win)/opt.n_shift+1);
      pup = nan(nseg,1);
      pup_df = nan(nseg,1);
      for j=1:nseg
        tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        pup(j) = mean(tmp.*gausswin(opt.n_win,3));
        if j==nseg
            tmp2 = pupil_df((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win-1);
            pup_df(j) = mean(tmp2.*gausswin(opt.n_win-1,3));
        else
            tmp2 = pupil_df((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
            pup_df(j) = mean(tmp2.*gausswin(opt.n_win,3));
        end
      end

      % -------------------------------
      % beamforming
      % -------------------------------
      para      = [];
      para.iscs = 1;
      para.reg  = 0.05;
      tp_filt   = tp_beamformer(real(tp_csd),lf,para);
      % -------------------------------
      % identify artifactual segments
      idx = find(~isnan(dataf(1,:))'&~isnan(pup)&~isnan(pup_df));
      % -------------------------------
      % sensor-level pupil power correlation
      % -------------------------------
      outp.sens_r(outp.chanidx>0,ifreq) = corr(pup(idx),abs(dataf(outp.chanidx(outp.chanidx>0),idx).^2)','type','spearman');
      outp.sens_r_df(outp.chanidx>0,ifreq) = corr(pup_df(idx),abs(dataf(outp.chanidx(outp.chanidx>0),idx).^2)','type','spearman');
      % -------------------------------
      % correlate power with pupil
      % -------------------------------
      src = abs(tp_filt'*dataf).^2; % source level power fluct
      outp.src_r(:,ifreq) = corr(pup(idx),src(:,idx)','type','spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx),src(:,idx)','type','spearman');
 
      clear src pup csd tp_csd dataf 
      
    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

error('!')
%%

for isubj = 1 : 37

load(sprintf('/home/tpfeffer/pp/proc/src/pp_mue_src_pupil_power_correlations_s%d_b1_v3.mat',isubj))

all(:,:,isubj) =outp.tp_src_r;

end

for igrid = 1 : max(BNA.tissue_5mm(:))
  all_BNA(igrid,:,:) = tanh(mean(atanh(all(BNA.tissue_5mm == igrid,:,:))));
end
