

ord    = pconn_randomization;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 4

for isubj = 1  : length(SUBJLIST)
for m = 1 : 3
  for iblock = 1 : 2
 

  im = find(ord(isubj,:)==m);

  load(sprintf('~/pconn/proc/preproc/pconn_rejected_comps_s%d_m%d_b%d_f1_v%d.mat',SUBJLIST(isubj),im,iblock,v))

  all_rej_comp(isubj,m,iblock) = sum(rej_comp);
  
  end
end
end