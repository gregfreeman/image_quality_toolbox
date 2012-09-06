function foveateWaveletSensitivityTest

w1=foveateWaveletSensitivity([512,512],[],4,1);
w3=foveateWaveletSensitivity([512,512],[],4,3);
w7=foveateWaveletSensitivity([512,512],[],4,7);
figure(1),imagesc(w3)
colormap gray
axis off
axis square
% title('\nu=3, p_f=[256 256]')
outputFigure(1,'wavelet_foveation_sensitivity',1,'-depsc2')

t=median(w1(:));
t
t=max(w1(:))*0.01;
t
m1=w1>t;

t=median(w3(:));
t
t=max(w3(:))*0.01;
t
m3=w3>t;

t=median(w7(:));
t
t=max(w7(:))*0.01;
t
m7=w7>t;

a1=sum(m1(:))/512^2;
a3=sum(m3(:))/512^2;
a7=sum(m7(:))/512^2;
a1
a3
a7



% figure(2),imagesc(m)
% colormap gray


figure(2)
colormap gray
subplot(3,1,1)
imagesc(w1)
title('v=1')
subplot(3,1,2)
imagesc(w3)
title('v=3')
subplot(3,1,3)
imagesc(w7)
title('v=7')


figure(3)
colormap gray
subplot(3,1,1)
imagesc(m1)
title('\nu=1')
subplot(3,1,2)
imagesc(m3)
title('v=3')
subplot(3,1,3)
imagesc(m7)
title('v=7')

sum(m1(:))
sum(m3(:))
sum(m7(:))

foveateWaveletSensitivity([512,512],[],4,7,struct('use_band_weighting',0,'plot_band_weighting',0));
foveateWaveletSensitivity([512,512],[],4,7,struct('use_band_weighting',0,'plot_band_weighting',1));
foveateWaveletSensitivity([512,512],[],4,7,struct('use_band_weighting',1,'plot_band_weighting',0));
foveateWaveletSensitivity([512,512],[],4,7,struct('use_band_weighting',1,'plot_band_weighting',1));

w=foveateSensitivity([512,512],[],3);
imagesc(w)
colorbar
