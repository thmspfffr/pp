%% pp_pupil_preproc
% preprocessing after pupil time series is converted in
% pp_pupil_timeseries.m

% 04/2019
clear

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
fsample         = 1000;

for isubj = SUBJLIST
  for m = 1:3
    
    fprintf('Processing subject %d, session%d ...\n',isubj,m)
    
    d       = dir(sprintf('~/pp/proc/pup/pp_pupil_diameter_s%d_m%d_b*.mat',isubj,m));
    d_evts  = dir(sprintf('~/pp/proc/pup/pp_pupil_events_s%d_m%d_b*.mat',isubj,m));

   for iblock = 1 : length(d)
     
     block = str2num(d(iblock).name(end-4));
      
     % load diameter timeseries (raw)
     load([d(iblock).folder '/' d(iblock).name])
     first_timestamp = pupil(1,1);
     
     % load events
     load([d_evts(iblock).folder '/' d_evts(iblock).name])
     blinksmp = blinks - first_timestamp + 1; clear blinks
     saccsmp  = saccs- first_timestamp + 1; clear saccs
     
     if length(pupil) < 200000
       fsample  = 250;
       blinksmp = round(blinksmp./4);
       saccsmp  = round(saccsmp./4);
     else
       fsample = 1000;
     end
     
     tmp = pupil;

     [pupil(:,4),newblinksmp,~,dat] = blink_interpolate(pupil, blinksmp, fsample, 0); 
     saccs=tp_detect_microsaccades(pupil(:,2:3),1000,12,6);
    
     pupil(:,5) = tmp(:,4);
     pupil(:,4) = blink_regressout(pupil(:,4), fsample, blinksmp, saccsmp, 0, 1);
%      pupil(:,4) = blink_regressout(pupil(:,4), fsample, saccs, saccs, 1, 1)
     pupil(:,2) = dat.gazex; % return the interpolated x gaze time series
     pupil(:,3) = dat.gazey; % return the interpolated y gaze time series
     art.saccs_ts = saccsmp;
     art.blinks_ts = blinksmp;
     
     if length(pupil) < 200000
       pupil=resample(pupil,4,1);
     end
     
     save(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,m,block),'art','pupil')
     clear pupil
   end
  end
end

error('!')
%% PLOT RESULTS

% d = dir('/home/tpfeffer/pp/proc/pup/pp_pupil_diameter_cleaned*');
% for dd = 1 : length(d)
%   load([d(dd).folder '/' d(dd).name])
%   plot(pupil(:,4))
%   drawnow
% end

