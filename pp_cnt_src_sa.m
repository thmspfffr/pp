%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% pp_cnt_src_sa

clear

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
grid_size = 'BNA';
v_rawdat = 6;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
% addpath /home/gnolte/neuconn/OLD/matlab/rest/
addpath(genpath('/home/gnolte/meth'));
addpath ~/pconn/matlab
ft_defaults

outdir = '/home/tpfeffer/pp/proc/src/';
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
%%
for im = 1 : 3
  for isubj = SUBJLIST
        %
    if ~exist(sprintf([outdir 'pp_cnt_sa_s%d_m%d_v%d_processing.txt'],isubj,im,v))
      system(['touch ' outdir sprintf('pp_cnt_sa_s%d_m%d_v%d_processing.txt',isubj,im,v)]);
    else
      continue
    end
    %
    fprintf('Processing s%d m%d ...\n',isubj,im);
    
    % ----- ---------------------------------------------------
    % SELECT MRI DATA
    % --------------------------------------------------------
    
    load sa_meg_template;
    
    fprintf('Looking for MRI ...\n');
    
    mridir = sprintf('/home/tpfeffer/pconn/rawdata/mri/p%d/mri/',isubj);
    c_mrid = dir(mridir);
    if isempty(c_mrid)
      mridir = sprintf('/home/tpfeffer/pconn/rawdata/mri/p%d/MRI/',isubj);
      c_mrid = dir(mridir);
    end
    
    br = 0;
        
    if ~isempty(c_mrid)
      if length(c_mrid) > 2
        for imri = 3 : length(c_mrid)
          if strcmp(c_mrid(imri).name(end-5:end),'V2.mri')
            mri_data = [mridir c_mrid(imri).name];
            fprintf('Looking for MRI ... Found!\n');
            br = 0;
          end
        end
      else
        fprintf('Looking for MRI ... Not Found!\n');
        br = 1;
      end
    else
      fprintf('Looking for MRI ... Not Found!\n');
      br = 1;
    end
    
    if br == 0
      
      aal = tp_aalgrid();
      sa_meg_template.grid_aal4mm = aal.grid_4mm;
      sa_meg_template.grid_aal6mm = aal.grid_6mm;

      vtpm = tp_create_grid('vtpm');
      sa_meg_template.grid_vtpm_4mm = vtpm.grid_4mm;
      sa_meg_template.grid_vtpm_6mm = vtpm.grid_6mm;
      
      if ~exist(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
        error('Create BNA file first')
      else
        load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
      end
      
      sa_meg_template.grid_BNA_5mm = BNA.grid_5mm;
      
      load ~/Documents/MATLAB/atlas_MMP1.0_4k.mat
      sa_meg_template.grid_MM = atlas.pos;
      
%       load ~/pmod/matlab/gene_values.mat
%       orig = [46 64 37];
      % Transform gene maps to MNI
%       locs=2*(locs-repmat(orig,[size(locs,1) 1]))/10;
%       sa_meg_template.grid_genemaps = locs;
      
%       idx = tp_aalgenemask(90);
%       sa_meg_template.grid_genemaps_aal = locs(idx,:);
      
      load cortex_brainstorm
      sa_meg_template.grid_cortex_brainstorm = cortex_brainstorm.pos;
      sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri_data);
      
    end
    
    saa = 0;
    
    % --------------------------------------------------------
    % SELECT SUBJECT & BLOCK
    % --------------------------------------------------------
    
    for iblock = 1 : 2
      
      fprintf('Processing block %d MEG-Data ...\n',iblock);
      fprintf('Loading MEG-Data ...\n');
      
      if br == 1
        %         mri_data = '/home/tpfeffer/pconn/rawdata/mri/p10/mri/aheitmann_V2.mri';
        %           disp(sprintf('Using template MRI ...'))
        load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
        load /home/gnolte/meth/templates/sa_template.mat
        
        aal = tp_aalgrid();
        sa_template.grid_aal4mm = aal.grid_4mm;
        sa_template.grid_aal6mm = aal.grid_6mm;
        
        %         m758 = tp_create_grid('m758');
        %         sa_template.grid_m758_4mm = m758.grid_4mm;
        %         sa_template.grid_m758_6mm = m758.grid_6mm;
        
        vtpm = tp_create_grid('vtpm');
        sa_template.grid_vtpm_4mm = vtpm.grid_4mm;
        sa_template.grid_vtpm_6mm = vtpm.grid_6mm;
        
        load ~/Documents/MATLAB/atlas_MMP1.0_4k.mat
        sa_template.grid_MM = atlas.pos;
        
        if ~exist(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
          error('Create BNA file first')
        else
          load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
        end
        
        sa_template.grid_BNA_5mm = BNA.grid_5mm;
        
%         load ~/pmod/matlab/gene_values.mat
%         orig = [46 64 37];
%         % Transform gene maps to MNI
%         locs=2*(locs-repmat(orig,[size(locs,1) 1]))/10;
%         sa_template.grid_genemaps = locs;
        
%         idx = tp_aalgenemask(90);
%         sa_template.grid_genemaps_aal = locs(idx,:);
        
        load cortex_brainstorm
        sa_template.grid_cortex_brainstorm = cortex_brainstorm.pos;
        
        sa  = tp_mk_sa_meg_withoutmri(sa_template,cfg1.headerfile);
        saa = 1;
        
      end
      
      clear meg_data cfg1 data cfg2
      
      load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      
      clear data
      meg_data = cfg1.headerfile;
      
      if ~saa
        sa 	= mk_sa_meg_forward(sa_meg1, meg_data);
      end
      
      if      strcmp(grid_size,'medium')
        
        L            = grid2L(sa.grid_medium_indi,sa.fp_indi);
        sa.L_medium  = L;
        sa.leadfield = 'medium';
        
      elseif  strcmp(grid_size,'coarse')
        
        L            = grid2L(sa.grid_coarse_indi,sa.fp_indi);
        sa.L_coarse  = L;
        sa.leadfield = 'coarse';
        
      elseif strcmp(grid_size,'cortex')
        
        L             = grid2L(sa.grid_cortex3000_indi,sa.fp_indi);
        sa.L_coarse   = L;
        sa.leadfield  = 'cortex';
        
      elseif strcmp(grid_size,'cortex_lowres')
        
        L             = grid2L(sa.grid_cortex_lowres_indi,sa.fp_indi);
        sa.L_coarse   = L;
        sa.leadfield  = 'cortex_lowres';
        
      elseif strcmp(grid_size,'cortex800')
        
        L             = grid2L(sa.grid_cortex800_indi,sa.fp_indi);
        sa.L_coarse   = L;
        sa.leadfield  = 'cortex800';
        
      elseif strcmp(grid_size,'xcoarse')
        
        L             = grid2L(sa.grid_xcoarse_indi,sa.fp_indi);
        sa.L_xcoarse  = L;
        sa.leadfield  = 'xcoarse';
        
      elseif strcmp(grid_size,'aal')
        
        load aalmask_grid_medium
        load sa_meg_template;
        para.grid = 'medium';
        pos = tp_grid2aal(sa.grid_medium_indi',para);
        sa.grid_aal_indi = pos';
        
        L             = grid2L(sa.grid_aal_indi,sa.fp_indi);
        sa.L_aal      = L;
        sa.leadfield  = 'aal';
        
      elseif strcmp(grid_size,'aal_6mm')
        
        L             = grid2L(sa.grid_aal6mm_indi,sa.fp_indi);
        sa.L_aal_6mm  = L;
        sa.leadfield  = 'aal_6mm';
        sa.aal_centr  = aal.centroids(1:91,:);
        sa.aal_label  = aal.tissue_6mm(aal.tissue_6mm>0 & aal.tissue_6mm<92);
        
      elseif strcmp(grid_size,'aal_4mm')
        
        L             = grid2L(sa.grid_aal4mm_indi,sa.fp_indi);
        sa.L_aal_4mm  = L;
        sa.leadfield  = 'aal_4mm';
        sa.aal_centr  = aal.centroids(1:91,:);
        sa.aal_label  = aal.tissue_4mm(aal.tissue_4mm>0 & aal.tissue_4mm<92);
        
      elseif strcmp(grid_size,'m758_4mm')
        
        L              = grid2L(sa.grid_m758_4mm_indi,sa.fp_indi);
        sa.L_m758_4mm  = L;
        sa.leadfield   = 'm758_4mm';
        
      elseif strcmp(grid_size,'m758_6mm')
        
        L              = grid2L(sa.grid_m758_6mm_indi,sa.fp_indi);
        sa.L_m758_6mm  = L;
        sa.leadfield   = 'm758_6mm';
        
      elseif strcmp(grid_size,'vtpm_4mm')
        
        L              = grid2L(sa.grid_vtpm_4mm_indi,sa.fp_indi);
        sa.L_vtpm_4mm  = L;
        sa.leadfield   = 'vtpm_4mm';
        
      elseif strcmp(grid_size,'vtpm_6mm')
        
        L              = grid2L(sa.grid_vtpm_6mm_indi,sa.fp_indi);
        sa.L_vtpm_6mm  = L;
        sa.leadfield   = 'vtpm_6mm';
        
      elseif strcmp(grid_size,'cortex_brainstorm')
        %
        L                       = grid2L(sa.grid_cortex_brainstorm_indi,sa.fp_indi);
        sa.L_cortex_brainstorm  = L;
        sa.leadfield            = 'cortex_brainstorm';
        
      elseif strcmp(grid_size,'genemaps')
        
        L = grid2L(sa.grid_genemaps_indi,sa.fp_indi);
        sa.L_genemaps           = L;
        sa.leadfield            = 'genemaps';
        
        
      elseif strcmp(grid_size,'genemaps_aal')
        
        L = grid2L(sa.grid_genemaps_aal_indi,sa.fp_indi);
        sa.L_genemaps_aal           = L;
        sa.leadfield            = 'genemaps_aal';
        
      elseif strcmp(grid_size,'MM')
        
        L = grid2L(sa.grid_MM_indi,sa.fp_indi);
        sa.L_genemaps_aal      	= L;
        sa.leadfield            = 'MM';
        
      elseif strcmp(grid_size,'BNA')
        
        L = grid2L(sa.grid_BNA_5mm_indi,sa.fp_indi);
        sa.L_genemaps_aal      	= L;
        sa.leadfield            = 'BNA_5mm';
        
      end
      
      save([outdir sprintf('pp_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v)],'sa');
      close all
    end
    br = 0;
  end
  
end

