x = pearsrnd(0,1,-0.9,4,1000,1);

y = -x + (rand(1000,1)-0.5)/10;

thresh = quantile(x,[0:0.05:1]);
idx = discretize(x,thresh);
   
    for i = 1 : 20
      mean_sig(i) = mean(y(idx==i));
    end
    
    
    thresh = linspace(min(x),max(x),20);
idx = discretize(x,thresh);
   
    for i = 1 : 20
      mean_sig2(i) = mean(y(idx==i));
    end
%       if isnan(mean(mean(mean(tmp_pxx(:,:,idx==i),3))))
%         a=1;

figure_w
subplot(2,3,1);

plot(x,y,'.')
axis square; tp_editplots
xlabel('x'); ylabel('y')
subplot(2,3,2);

plot(mean_sig)
axis square; tp_editplots
xlabel('Pupil bin'); ylabel('y')

subplot(2,3,3);

plot(mean_sig2)
axis square; tp_editplots
xlabel('Pupil bin'); ylabel('y')

y = -x.^2 - x + (rand(1000,1)-0.5)/10;

thresh = quantile(x,[0:0.05:1]);
idx = discretize(x,thresh);
   
    for i = 1 : 20
      mean_sig(i) = mean(y(idx==i));
    end
    
    
    thresh = linspace(min(x),max(x),20);
idx = discretize(x,thresh);
   
    for i = 1 : 20
      mean_sig2(i) = mean(y(idx==i));
    end
%       if isnan(mean(mean(mean(tmp_pxx(:,:,idx==i),3))))
%         a=1;

% figure_w
subplot(2,3,4);

plot(x,y,'.')
axis square; tp_editplots
xlabel('x'); ylabel('y')
subplot(2,3,5);

plot(mean_sig)
axis square; tp_editplots
xlabel('Pupil bin'); ylabel('y')

subplot(2,3,6);

plot(mean_sig2)
axis square; tp_editplots
xlabel('Pupil bin'); ylabel('y')

