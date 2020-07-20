%% pp_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
ord = pconn_randomization;
v = 1;

addpath ~/pconn/matlab
for isubj = SUBJLIST
  
  fprintf('Processing subject %d...\n',isubj)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    
    % load meg data
    load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
    % load pupil data
    load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))

    % filter pupil time series
    % bandpass filter: 0.005 - 2 hz
    % gets rid of slow drift
    k = 2;
    f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil = resample(pupil,400,1000);
    
    data.trial = data.trial(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    len = min([size(pupil,1) size(data.trial,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    data.trial = data.trial(:,1:len);
    pupil = pupil(1:len);
    
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    % replace this later
    pupil(isnan(data.trial(1,:)))=[];
    data.trial(:,isnan(data.trial(1,:))) = [];
    
    flp = 8;           % lowpass frequency of filter
    fhi = 12;
    
    para.ord = 4;
    k=4;                  % 2nd order butterworth filter
    fnq=400/2;       % Nyquist frequency
    Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    [bfilt,afilt]=butter(k,Wn);
    
    data_f = abs(hilbert(single(filtfilt(bfilt,afilt,data.trial'))));
    
    k = 2;
    f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    data_f = filtfilt(bhil, ahil, double(data_f));

    
    if isubj >= 32
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
    elseif isubj < 4
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
    elseif isubj == 17
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
    else
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,im))
    end
          
%     load /home/tpfeffer/pconn/proc/src/pconn_sa_s26_m2_b1_v4.mat
%     pars            = [];
%     pars.markersize = 0;
%     pars.linewidth  = 9;
%     pars.cbar       = 0;
%     pars.scale      = [-0.05 0.05]
%     pars.cmap       = jet;
% 
%     figure; set (gcf,'color','w')
    [tmp_r ,p]=corr(data_f,pupil);
    r(:,isubj,iblock) = pconn_sens_interp274(idx,tmp_r);
%     pars.resolution = 600;
%     showfield_colormap(r(:,isubj,iblock),sa.locs_2D,pars);
%     drawnow
    
    
    
    
  end
end

