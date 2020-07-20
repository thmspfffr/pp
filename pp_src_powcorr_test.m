
%% pp_src_powcorr_test

% Call pupmod_all_powcorr_periphereal.m next, in order to clean
% estimated FCs from artifacts.

clear
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v               = 23;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'cortex_lowres';
% foi_range       = 2.^[1:.25:7];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% fsample         = 400;
% segleng         = 80;
% segshift        = 40;
% width           = 4;
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v               = 24;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% foi_range       = 2.^[1:.25:7];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% weighting       = 1;
% width           = 4;
% --------------------------------------------------------
% VERSION 25 - VOXEL LEVEL, 90 samples AAL (6mm sampling)
% --------------------------------------------------------
% v               = 25;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_6mm';
% foi_range       = 2.^[1:.25:6];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% weighting       = 0;
% width           = 4;
% FOI             = [10 11 12 13];
% --------------------------------------------------------
% VERSION 25 - VOXEL LEVEL, 90 samples AAL (6mm sampling)
% --------------------------------------------------------
% v               = 35;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_6mm';
% foi_range       = 2.^[1:.25:6];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% weighting       = 0;
% width           = 4;
% FOI             = [10];
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v               = 33;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'cortex_lowres';
% foi_range       = 2.^[1:.25:7];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% fsample         = 400;
% segleng         = 80;
% segshift        = 40;
% width           = 4;
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v               = 43;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
foi_range       = 2.^[1:.25:7];
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 0;
allpara.tau     = nan;
fsample         = 400;
segleng         = 80;
segshift        = 40;
width           = 4;
% --------------------------------------------------------

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid, 'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'genemaps')
  v_grid = 13;
elseif strcmp(allpara.grid,'genemaps_aal')
  v_grid = 14;
elseif strcmp(allpara.grid,'cortex800')
  v_grid = 16;
end

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/
ft_defaults()
outdir   = '/home/tpfeffer/pp/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/

ord   = pconn_randomization;
%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
%     
%     if ~exist(sprintf([outdir 'pp_src_powcorr_test_s%d_m%d_v%d_processing.txt'],isubj,m,v))
%       system(['touch ' outdir sprintf('pp_src_powcorr_test_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%     else
%       continue
%     end

%     if ~weighting
      if strcmp(allpara.grid,'cortex_lowres')
        powcorr = zeros(400,400,2,length(foi_range));
      else
        powcorr = zeros(90,90,2,length(foi_range));
      end
%     else
%       powcorr = zeros(90,90,2,length(foi_range));
%     end
    for iblock = 1:2
      
      pars = [];
      pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      sa        = load(pars.sa);
      
      fprintf('Loading MEG data ...\n');
      
      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        if ~exist(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
          if strcmp(allpara.grid,'cortex_lowres')
            powcorr(:,:,iblock,:) = nan(400,400,1,length(foi_range));
          else
            powcorr(:,:,iblock,:) = nan(90,90,1,length(foi_range));
          end
          continue
        else
          error('Data corrupt?')
        end
      end
      %       data.trial{1} = data.trial;
      %       dat = megdata2mydata(data);
      if isempty(data.start_of_recording) && ~isempty(data.end_of_recording)
        if (data.end_of_recording-600*data.fsample)<1
          data.start_of_recording = 1;
        else
          data.start_of_recording = data.end_of_recording-600*data.fsample;
        end
      elseif ~isempty(data.start_of_recording) && isempty(data.end_of_recording)
        if (data.start_of_recording+600*data.fsample)>size(data.trial,2)
          data.end_of_recording = size(data.trial,2);
        else
          data.end_of_recording = data.start_of_recording+600*data.fsample;
        end
      elseif isempty(data.start_of_recording) && isempty(data.end_of_recording)
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
      elseif isnan(data.start_of_recording) && isnan(data.end_of_recording)
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
      end
      
      data.trial = data.trial(:,data.start_of_recording:data.end_of_recording);
      data.time = data.time(data.start_of_recording:data.end_of_recording);
      
      for ifoi = 1:length(foi_range)
        
        fprintf('Processing s%d m%d b%d f%d  ...\n', isubj,m,iblock,ifoi)
%         data.time(isnan(data.trial(1,:)))=[];
%         data.trial(:,isnan(data.trial(1,:)))=[];
        
        para = [];
        para.freq = foi_range(ifoi);
        para.fsample = 400;
%         cs = tp_wavelet_crossspec(data,foi_range(ifoi),width);
        cs = tp_compute_csd_wavelets(data,para);
        
        para.iscs = 1;
        para.reg  = 0.05;
        if strcmp(allpara.grid,'cortex_lowres')
          sa.sa.filt      = pconn_beamformer(real(cs),sa.sa.L_coarse,para);
        else
          sa.sa.filt      = pconn_beamformer(real(cs),sa.sa.L_aal_6mm,para);
        end
        fprintf('Computing pow corr  ...\n')
          
        powcorr(:,:,iblock,ifoi) = tp_data2orthopowcorr_wavelet(data,foi_range(ifoi),sa);       
        
      end
           
    end
    save(sprintf([outdir 'pp_src_powcorr_test_s%d_m%d_v%d.mat'],isubj,m,v),'powcorr');
    clear powcorr

  end
end
error('!')

%











%%
ord = pconn_randomization;

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    
    load(sprintf('/home/tpfeffer/pp/proc/conn/pp_cnt_src_powcorr_test_s%d_m%d_v43.mat',isubj,im))
    
    allpow(:,:,isubj,m,2,:) = squeeze(nanmean(powcorr,3));
    
    load(sprintf('/home/tpfeffer/pp/proc/conn/pp_src_powcorr_test_s%d_m%d_v43.mat',isubj,im))
    
    allpow(:,:,isubj,m,1,:) = squeeze(nanmean(powcorr,3));
    
  end
end
    

allpow = allpow(:,:,SUBJLIST,:,:,:);

