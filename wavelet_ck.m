function meg_preproc_step2_ck(sourcedir,subjCode,subjMegrun,gridRes)

%%% time-frequency decomposition & source projection

% where to find the data:
fileprfx=[sourcedir,subjCode,subjMegrun];

% 1.1 read in the preprocessed and combined MEG & pupil data
load([fileprfx 'meg_preproc.mat'],'meg');
data=meg; meg=[];

% 1.2 read in original file header to retrieve sensor info
sensloc=ft_read_header([fileprfx 'c,rfDC']);

%%% 1.3 loop over frequencies & source reconstruct
freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci

% recover sampling rate from time axis
srate    =round(1/diff(data.time{1}([1 2]),[],2));

% highpass filter on MEG (not! eye data)
fprintf('Applying highpass filter to MEG data...\n')
megChan=~cellfun(@isempty,regexp(data.label,'A*'));
data.trial{1}(:,megChan)=ft_preproc_highpassfilter(data.trial{1}(:,megChan),srate,1,4,'but');

% % select pupil data
% eye=ft_selectdata(struct('channel',{'Eye*'}),data);
% % compute 1st order derivative as an extra channel
% eye.trial{1}(end+1,:)=gradient(eye.trial{1}(strcmp(eye.label,'EyePupilFilt'),:));
% eye.label{end+1}='EyePupilFiltGrad';

%% spectral decomp and spatial filter computation
for ifreq=1:numel(freqoi)
% freq analysis, pass 1, yields TF Fourier rep ######################
% this is for multiplication with the ortho-normalised spatial filter
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
% reset freq to requested freq
tf.freq=freqoi(ifreq);

% mark time bins that coincide with relevant artifacts - outdated
% artifPnts=sortrows([data.cfg.artfctdef.visual.artifact;
%                     data.cfg.artfctdef.blink.artifact]);

% blinks have been removed via ICA, this is mostly myogenic, sq jumps
artifPnts=data.cfg.artfctdef.visual.artifact;
tfPnts   =tf.time*srate;
% only discard those bins that "center" on an artifact
critDist =diff(tfPnts([1,2]),[],2)/2; % set critical dist to central 50%
keepBins=logical([]);
% discard TF bins that overlap with artifacts (set zero in keepBins)
for ibin=1:numel(tfPnts)
    keepBins(ibin)=~any(abs(tfPnts(ibin)-ceil(mean(artifPnts,2)))<critDist);   
end    

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
csdData=tf.fourierspctrm(:,:,:,keepBins);
for itbin=1:numel(csdTime)
    fspec=squeeze(csdData(:,:,:,itbin)).';
    for ichan=1:numel(tf.label);
        csd(:,ichan)=csd(:,ichan)+fspec(ichan)*conj(fspec);
    end
end
csd=csd./numel(tf.time); % avg cross-spectral dens matrix
csdData=[]; csdTime=[];
% #### end of pass 1

% % subsample pupil data according to tf temp resolution
% resEye=ft_resampledata(struct('time',{{tf.time}}),eye);

% freq analysis, pass 2, yields TF pow and CSD rep ######################
% creates input data for fieldtrip DICS implementation
cfg.output='powandcsd';
tfpow=ft_freqanalysis(cfg,data);
tfpow.freq=freqoi(ifreq);

% exclude artifacts from filter computation
tfpow.powspctrm=tfpow.powspctrm(:,:,keepBins);
tfpow.crsspctrm=tfpow.crsspctrm(:,:,keepBins);
tfpow.time=tfpow.time(1:sum(keepBins));

% load indiv leadfield and headmodel
srcmod_grid=load(sprintf('%s%s/grid%d',sourcedir,subjCode,gridRes));
srcmod_vol =load(sprintf('%s%s/vol%d',sourcedir,subjCode,gridRes));

%%% check on alignment - looking good
% figure
% ft_plot_mesh(srcmod_grid.indgrid.pos(srcmod_grid.indgrid.inside,:)), hold on
% ft_plot_vol(ft_convert_units(srcmod_vol.vol,'m'),'facealpha',.5);
% ft_plot_sens(sensloc.grad,'markercolor','r','edgealpha',.5)

cfg=[];
cfg.method='dics';
cfg.frequency=freqoi(ifreq);
cfg.grid=srcmod_grid.indgrid;
cfg.headmodel=ft_convert_units(srcmod_vol.vol,'m');
cfg.grad = sensloc.grad;
cfg.channel={'MEG'};
cfg.dics.lambda='5%'; % standard value, think about individual optim
cfg.dics.keepfilter='yes';
cfg.dics.feedback='none';
cfg.dics.realfilter='yes';
tmp=ft_sourceanalysis(cfg,tfpow);

% grab spatial filter from FT output
realFilt=cat(3,tmp.avg.filter{:});
xF=nan(size(realFilt,2),size(realFilt,3));

% linearly combine orthogonal filters & normalise leadfield
for idip=1:size(realFilt,3) % loop over voxels/dipoles
    lF=realFilt(:,:,idip);
    % normalize leadfields
    %        lFnorm=norm(lF.', 'fro'); FT version 1
    %        lFnorm=sum(lF(:).^2)^.5;  FT version 2 - .5 = 'normparam'
    lFnorm=sum(lF(:).^2); % as per Hipp et al. Nat Neu 2012
    lF=lF./lFnorm;
    fW=eig(lF*real(csd)*lF.');
    % lin combination
    xF(:,idip)=fW(1)*squeeze(realFilt(1,:,idip)).'+...
               fW(2)*squeeze(realFilt(2,:,idip)).'+...
               fW(3)*squeeze(realFilt(3,:,idip)).';
end

%%% non-essential => has been moved to next step to reduce storage req
% pass data through spatial filter:
%beamedTF=xF.'*squeeze(tf.fourierspctrm(1,:,1,:));
%realFilt=[];

source=[];
source.tfDat=squeeze(tf.fourierspctrm(1,:,1,:));
source.spatFilt=xF.';
source.time=tf.time;
source.cleanBin=keepBins; % store artifact info for later use
source.freq=freqoi(ifreq);
% label voxels
source.label=cellfun(@(x) sprintf('vxl%d',x),num2cell(1:sum(tmp.inside)).','UniformOutput',false);
source.dimord='non_standard';
tf=[]; tfpow=[]; tmp=[]; xF=[];

% sourceTF.trial(:,1,:)=permute(beamedTF,[1,3,2]);
% % calculate corr between orthogon pow envelopes (Hipp et al 2012 NN)
% sourceTF.corrOrtho=ft_connectivity_powcorr_ortho(beamedTF);
% beamedTF=[];

% % add eye chan label
% addchan=1:numel(eye.label); % aux var
% source.label(end+addchan)=eye.label;
% %source.trial(end+addchan,1,:)=cat(3,eye.trial{:});
% source.eyeDat=resEye.trial{:};

% % add source power envelopes
% sourceTF.powEnv=20*log10(abs(sourceTF.trial));
% sourceTF.corrPupil=corr(sourceTF.powEnv.',eye.trial{1}(4,:).');

xfreq{ifreq}=source;
end
save([fileprfx 'meg_specest_canonF_wavelet_source_' num2str(gridRes) 'mm.mat'],'xfreq','-v7.3');
return