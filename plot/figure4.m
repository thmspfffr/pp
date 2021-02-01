%% PLOT SOURCE MAPS: TASK VS REST
% load /home/gnolte/meth/templates/mri.mat
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

[plt_gla,plt_hh,plt_mue,plt_hh_cnt,plt_all]=pp_load_results(2);

for ifoi = [5 11 14 22]
  
  figure_w
  
  [h,p] = ttest(plt_all.corr_src(:,ifoi,:),zeros(size(plt_all.corr_src(:,ifoi,:))),'dim',3);
  h=p<fdr1(p(:),0.05,1);
  par=nanmean(plt_all.corr_src(:,ifoi,:),3).*h;

  % project onto fine grid
  par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);
  
  
  clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
  para = [];
  para.colorlimits = clim
  para.colormaps{1} = cmap;
  para.orientation = 'axial';
  
  para.dslice_shown = 0.95;
  para.colorbar= 0;
  
  tp_showmri_transp(mri,para,[sa_template.grid_fine par]); f=get(gcf)
  
  move_x = [0:0.08:0.32];
  move_y = [0.08:-0.02:0];
  
  for i=4:-1:1
    for j = 0:3
      if j ~= 4
        f.Children(i*4-j).Position(1)=f.Children(i*4-j).Position(1)-move_x(j+1);
      end
      f.Children(i*4-j).Position(2)=f.Children(i*4-j).Position(2)+move_y(i);
    end
  end
  
  text(1,1,sprintf('[%.3f %.3f]\n [%.3f Hz]',clim(1),clim(2),ifoi))  

  set(gcf,'renderer','painters')
  print(gcf,'-dpdf',sprintf('~/pp/plots/pp_figure4_dt0_f%d_v%d.tiff',ifoi,v))
  
end

[plt_gla,plt_hh,plt_mue,plt_hh_cnt,plt_all]=pp_load_results(1);

for ifoi = [5 11 14 22]
  
  figure_w
  
  [h,p] = ttest(plt_all.corr_src_df(:,ifoi,:),zeros(size(plt_all.corr_src_df(:,ifoi,:))),'dim',3);
  h=p<fdr1(p(:),0.05,1);
  par=nanmean(plt_all.corr_src_df(:,ifoi,:),3).*h;

  % project onto fine grid
  par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);
  
  
  clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
  para = [];
  para.colorlimits = clim
  para.colormaps{1} = cmap;
  para.orientation = 'axial';
  
  para.dslice_shown = 0.95;
  para.colorbar= 0;
  
  tp_showmri_transp(mri,para,[sa_template.grid_fine par]); f=get(gcf)
  
  move_x = [0:0.08:0.32];
  move_y = [0.08:-0.02:0];
  
  for i=4:-1:1
    for j = 0:3
      if j ~= 4
        f.Children(i*4-j).Position(1)=f.Children(i*4-j).Position(1)-move_x(j+1);
      end
      f.Children(i*4-j).Position(2)=f.Children(i*4-j).Position(2)+move_y(i);
    end
  end
  
  text(1,1,sprintf('[%.3f %.3f]\n [%.3f Hz]',clim(1),clim(2),ifoi))  

  set(gcf,'renderer','painters')
  print(gcf,'-dpdf',sprintf('~/pp/plots/pp_figure4_dt1_f%d_v%d.tiff',ifoi,v))
  
end