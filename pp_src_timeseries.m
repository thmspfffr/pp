%% pp_src_timeseries
% compute amplitude/power time series in source space
% source models are loaded from pp_create_sourcemodel.m
% wavelets based on christian keitels code
% 04/2019

v       = 1;
freqoi  = 2.^(1:(1/4):7);

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults


for isubj = SUBJLIST
  for m = 1 : 3
    for iblock = 1 : 2
      
      % insfert alenas data here
      load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      
      % 1.2 read in original file header to retrieve sensor info
%       sensloc=ft_read_header([fileprfx 'c,rfDC']);
      
      %%% 1.3 loop over frequencies & source reconstruct
       % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci
      srate    = 400;
      for ifreq=1:numel(freqoi)
        
        cfg=[];
        cfg.method='wavelet';
        cfg.output='fourier';
        cfg.channel={'MEG'};
        cfg.foi=freqoi(ifreq);
        cfg.width=5.83; % again, as per Hipp et al. (2012) Nat Neurosci
        tempSD=1./(2*pi*(cfg.foi./cfg.width)); % temporal SD in sec
        tempWd=round(3*tempSD*srate)/srate; % set to 3 std dev, comp with 1000 Hz sampl rate
        cfg.toi=tempWd.*(1:floor(data.time{1}(end)./tempWd));
        
        % keep in mind that fieldtrip uses a proprietary setting of the gwidth
        % parameter (default = 3 cycles) internally that is independent of the
        % here requested time axis
        cfg.pad='nextpow2';
        tf=ft_freqanalysis(cfg,data);
        
 
        % ARTIFACT REMOVAL
%         artifPnts=data.cfg.artfctdef.visual.artifact;
%         tfPnts   =tf.time*srate;
%         critDist =diff(tfPnts([1,2]),[],2)/2; % set critical dist to central 50%
%         keepBins=logical([]);
%         for ibin=1:numel(tfPnts)
%           keepBins(ibin)=~any(abs(tfPnts(ibin)-ceil(mean(artifPnts,2)))<critDist);
%         end
        
        % additionally omit edge datapoints to excl artif
        timeAx=tf.time; % original time axis
        timeKp=dsearchn(timeAx.',[3;timeAx(end)-3]); % excl ~ 1st & last 3 sec
        keepBins(1:timeKp(1)-1)=false;
        keepBins(timeKp(2):end)=false;
        
        % compute cross-specral density by time bin, then average csd's
        % ### this is a workaround to have the CSD for extra LF normalisation
        csd=zeros(numel(tf.label)*[1 1]);
        % csd calc excludes artifact bins
        csdTime=tf.time(keepBins);
        csdData=tf.fourierspctrm(:,:,:,1:end-1);
        for itbin=1:numel(csdTime)-1
          fspec=squeeze(csdData(:,:,:,itbin)).';
          for ichan=1:numel(tf.label);
            csd(:,ichan)=csd(:,ichan)+fspec(ichan)*conj(fspec);
          end
        end
        csd=csd./(numel(tf.time)-1); % avg cross-spectral dens matrix
        csdData=[]; csdTime=[];
        % #### end of pass 1
         
        %% SOURCE ANALYSIS
        % -----------------------------------
        
        load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v)],'sa');
        
        para      = [];
        v         = 1;
        para.iscs = 1;
        para.reg  = 0.05;
        filt      = tp_beamformer(csd,sa.L_BNA_5mm,para);
        
        % Project data into source space
        % -----------------------------------
        
      end
    end
  end
  
  
  
  
  
  
  
