%% pp_asc2mat_cnt
% converts pupil area into diameter timeseries
% extracts events from asc file
% extracts time series for counting and button conditions

% 04/2019

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% Load file which states which block is count and which is button
% script in ~/pp/matlab/*oragnize*
% 0 = rest; 1 = counting; 2 = pressing
load('~/pp/proc/blocks_taskrest.mat');

for isubj = SUBJLIST
  for m = 1:3
    
    fprintf('Processing subject %d, session%d ...\n',isubj,m)
    
    d = dir(sprintf('/home/tpfeffer/pconn/rawdata/edf/p%d/s%d/edf_samples/*asc',isubj,m));
    d_evts  = dir(sprintf('/home/tpfeffer/pconn/rawdata/edf/p%d/s%d/edf_events/*asc',isubj,m));
    
    % note that cond is only for subjects 1-28, i.e., those subjects that
    % were included in the final analysis.
    curr_cond = cond(:,isubj,m);
    
    cnt   = find(curr_cond==1);
    bttn  = find(curr_cond==2);
    
    if length(d)~=length(d_evts)
      error('!')
    end
    
    % FIRST EXTRACT COUNTING
    for icnt = 1 : length(cnt)
      if icnt == 1
        block = 1;
      else
        block = 2;
      end
      
      fprintf('Processing counting block%d...\n',block)

      for idir = 1 : length(d)     
        if ~isempty(strfind(d(idir).name,sprintf('b%d',cnt(icnt))))     
          break   
        end
      end
      
      if isempty(idir)
        continue
      end
            
      % ---------------
      % READ PUPIL DATA FROM *ASC FILE
      % ---------------
      fid = fopen([d(idir).folder '/' d(idir).name]);
      data=textscan(fid,'%q%q%q%q%*q','headerlines',0,'delimiter','\t');
      fclose(fid);
      
      fprintf('Reading data from ascii-file...\n')
      % Speed improvements possible here!
      for icell = 1 : numel(data)
        fprintf('Reading data from ascii-file... %d/%d...\n',icell,numel(data))
        if icell == 1
          C = data{icell};
          S = sprintf('%s*', C{:});
          samples(:,icell) = sscanf(S, '%f*')';
        else
          samples(:,icell)=zeros(length(data{icell}),1);
          samples(:,icell)=str2double(data{icell});
        end
      end
      clear pupil
      pupil  = samples;
      clear samples
      % DETERMINE WHETHER PUPIL AREA OR DIAMETER WAS RECORDED
      % type == 1: area | type == 2: diameter | type == 0: sth else
      % ---------------
      
      fn = fopen([d_evts(idir).folder '/' d_evts(idir).name]);
      while 1
        tmp = fgetl(fn);
        tmp = tmp(find(~isspace(tmp)));
        if ~isempty(findstr(tmp,'PUPILAREA'))
          type = 1;
          break
        elseif ~isempty(findstr(tmp,'PUPILDIAMETER'))
          type = 2;
          break
        end
        % if neither 1 nor 2, set 0
        type = 0;
      end
      fclose(fn);
      if type == 1
        fprintf('Determined type: AREA ...\n')
      else
        fprintf('Determined type: DIAMETER ...\n')
      end
      
      % ---------------
      % if pupil area, convert to diameter
      if type == 1
        pupil(:,4) = 256.*sqrt(pupil(:,4)./pi);
      end
      
      % save pupil time series      
      save(sprintf('~/pp/proc/pup/pp_pupil_diameter_counting_s%d_m%d_b%d.mat',isubj,m,block),'pupil')

      fid = fopen([d_evts(idir).folder '/' d_evts(idir).name]);

      % -------------------------------
      % EXTRACT SACCADES AND BLINKS
      % -------------------------------
      saccs = []; blinks = [];
      i=0;j=0;l=0;k=0;trlNum=0;
      while ~feof(fid)     
        tmp = fgetl(fid);
        if ~isempty(strfind(tmp,'pretrial_sync'))
          i = i + 1;
          tmp_trl(i,1)=cell2mat(textscan(tmp, '%*s%d%*[\n]'));
        elseif ~isempty(strfind(tmp,'posttrial_sync'))
          j = j + 1;
          tmp_trl(j,2)=cell2mat(textscan(tmp, '%*s%d%*[\n]'));
        elseif ~isempty(strfind(tmp, 'EBLINK'))
          k = k + 1;
          buff = textscan(tmp, '%s%s%d%d%*[\n]');
          blinks(k,1:2)=cell2mat(buff(3:4));
          clear buff
        elseif ~isempty(strfind(tmp, 'ESACC'))
          l = l + 1;
          buff = textscan(tmp, '%s%s%d%d%*[\n]');
          saccs(l,1:2)=cell2mat(buff(3:4));
          clear buff      
        end
      end
      % -------------------------------

      % save events 
    	save(sprintf('~/pp/proc/pup/pp_pupil_events_counting_s%d_m%d_b%d.mat',isubj,m,block),'saccs','blinks')
  
    end
    
    % NOW EXTRACT BUTTON
    for icnt = 1 : length(bttn)
      
      if isempty(idir)
        continue
      end
           
      if icnt == 1
        block = 1;
      else
        block = 2;
      end
      
      fprintf('Processing button block%d...\n',block)

      for idir = 1 : length(d)     
        if ~isempty(strfind(d(idir).name,sprintf('b%d',bttn(icnt))))     
          break   
        end
      end
            
      % ---------------
      % READ PUPIL DATA FROM *ASC FILE
      % ---------------
      fid = fopen([d(idir).folder '/' d(idir).name]);
      data=textscan(fid,'%q%q%q%q%*q','headerlines',0,'delimiter','\t');
      fclose(fid);
      
      fprintf('Reading data from ascii-file...\n')
      % Speed improvements possible here!
      for icell = 1 : numel(data)
        fprintf('Reading data from ascii-file... %d/%d...\n',icell,numel(data))
        if icell == 1
          C = data{icell};
          S = sprintf('%s*', C{:});
          samples(:,icell) = sscanf(S, '%f*')';
        else
          samples(:,icell)=zeros(length(data{icell}),1);
          samples(:,icell)=str2double(data{icell});
        end
      end
      clear pupil
      pupil  = samples;
      clear samples
      % DETERMINE WHETHER PUPIL AREA OR DIAMETER WAS RECORDED
      % type == 1: area | type == 2: diameter | type == 0: sth else
      % ---------------
      
      fn = fopen([d_evts(idir).folder '/' d_evts(idir).name]);
      while 1
        tmp = fgetl(fn);
        tmp = tmp(find(~isspace(tmp)));
        if ~isempty(findstr(tmp,'PUPILAREA'))
          type = 1;
          break
        elseif ~isempty(findstr(tmp,'PUPILDIAMETER'))
          type = 2;
          break
        end
        % if neither 1 nor 2, set 0
        type = 0;
      end
      fclose(fn);
      if type == 1
        fprintf('Determined type: AREA ...\n')
      else
        fprintf('Determined type: DIAMETER ...\n')
      end
      
      % ---------------
      % if pupil area, convert to diameter
      if type == 1
        pupil(:,4) = 256.*sqrt(pupil(:,4)./pi);
      end
      
      % save pupil time series      
      save(sprintf('~/pp/proc/pup/pp_pupil_diameter_button_s%d_m%d_b%d.mat',isubj,m,block),'pupil')

      fid = fopen([d_evts(idir).folder '/' d_evts(idir).name]);

      % -------------------------------
      % EXTRACT SACCADES AND BLINKS
      % -------------------------------
      saccs = []; blinks = [];
      i=0;j=0;l=0;k=0;trlNum=0;
      while ~feof(fid)     
        tmp = fgetl(fid);
        if ~isempty(strfind(tmp,'pretrial_sync'))
          i = i + 1;
          tmp_trl(i,1)=cell2mat(textscan(tmp, '%*s%d%*[\n]'));
        elseif ~isempty(strfind(tmp,'posttrial_sync'))
          j = j + 1;
          tmp_trl(j,2)=cell2mat(textscan(tmp, '%*s%d%*[\n]'));
        elseif ~isempty(strfind(tmp, 'EBLINK'))
          k = k + 1;
          buff = textscan(tmp, '%s%s%d%d%*[\n]');
          blinks(k,1:2)=cell2mat(buff(3:4));
          clear buff
        elseif ~isempty(strfind(tmp, 'ESACC'))
          l = l + 1;
          buff = textscan(tmp, '%s%s%d%d%*[\n]');
          saccs(l,1:2)=cell2mat(buff(3:4));
          clear buff      
        end
      end
      % -------------------------------

      % save events 
    	save(sprintf('~/pp/proc/pup/pp_pupil_events_button_s%d_m%d_b%d.mat',isubj,m,block),'saccs','blinks')
  
    end
    
  end
end







