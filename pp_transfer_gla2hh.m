function trans = pp_transfer_gla2hh()

load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
load ~/standard_sourcemodel_BNA_5mm.mat

pos_gla=sourcemodel.pos(sourcemodel.inside,:);
pos_hh=BNA.grid_5mm;

for i = 1 : 8799
  trans(i)=find(sum(pos_hh(i,:)==pos_gla,2)==3);
end

