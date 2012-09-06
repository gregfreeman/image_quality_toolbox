


%% test that the filter frequency response has the appropriate cutoff
N=10;
figure(1)
hold off
for i=1:N
    sigma=sqrt(2*log(2))/pi*N/i;
    sigmas(i)=sigma;
    l=ceil(sigma*3+1);
    t1=-l:l;
    h1=exp(-(t1.^2)./2./sigma./sigma)./sqrt(2*pi*sigma*sigma);
    sum(h1)
    h1=h1./sum(h1);
    % fvtool(h1,[1])
    n=512;
    h2=[h1 zeros(1,n-length(h1))];
    f=linspace(-1,1,n+1);
    f=f(1:n);
    plot(f,fftshift(abs(fft(h2)).^2))
    if i==1
        hold on
    end

    pass=fftshift(abs(fft(h2)).^2)>0.5;
    epass=fftshift(abs(fft(h2)))>0.5;
    fc=max(abs(f(pass)));
    fc2=max(abs(f(epass)));
    fcs(i)=fc;
    fcs2(i)=fc2;
end

fcs=fcs(1:end-2);
fcs2=fcs2(1:end-2);
sigmas=sigmas(1:end-2);

figure(2)
% hold off
plot(fcs,1./sigmas,'b')
xlabel 'freq (cutoff)'
ylabel '1/\sigma'


ratio=sqrt(2*log(2))/pi;

figure(3)
% hold off
plot(fcs.*sigmas.*pi./sqrt(2*log(2)),'b*-')
xlabel 'case'
ylabel 'ratio'

figure(4)
% hold off
plot(fcs2.*sigmas.*pi./sqrt(2*log(2)),'b*-')
xlabel 'case'
ylabel 'ratio'

% hold on
% plot(fcs,sigmas.^2,'g')
% plot(fcs.^2,sigmas,'r')
% plot(sqrt(fcs),sigmas,'c')
% plot(fcs,sqrt(sigmas),'k')

%% compare different models of foveation rolloff

x=(-n/2):(n/2-1);

platSize = n/4;
minBandwidth = 0.1;
HPBW = n/4;
maskBins = minBandwidth .*ones(size(x));
lambda = log(2)/HPBW;
r1=@(x) sqrt((x).^2);
r2=@(x) sqrt((x).^2)-platSize;


curr = r2(x);
curr(curr<0) =0;
fovControl = exp(-curr*lambda);

maskBins(fovControl>maskBins) = fovControl(fovControl>maskBins);

figure(5)
hold off
h1=plot(maskBins.*0.5,'g')

numBands=8;
%Use this number in areas where it exceeds the existing bandwidth
maskIndex=ceil(maskBins*numBands);
maskIndex(maskIndex<1)=1;
maskIndex(maskIndex>numBands)=numBands;

d = r1(x);

view_dist=1;
alpha=0.106;
e_2=2.3;
CT_0=1/64;

e=atan(d./view_dist./n).*180./pi; % eccentricity in degrees.

% f[i] = s/pow(sep,-i)*0.0175*viewDistance/pixelWidth;
% norm = max(norm,exp(-alpha*f[i]));
% e[y][x] = min(e[y][x],atan(sqrt(pow(fovea[i]*a.cols-x,2)+pow(fovea[numFovea+i]*a.rows-y,2))*pixelWidth/viewDistance)*180.0/PI);			
% double w = exp(-alpha*f[i]*(e[y+5][x+5]+e2)/e2)/norm;
% e_prime = atan(d*pixelWidth/viewDistance)*180.0/pi;			

% Sf_prime = exp(-alpha.*f.*(e+e_2)./e_2);

%CT=CT_0*exp(alpha.*f.*(e+e_2)./e_2);
%CS=1./CT;
f_c=e_2*log(1/CT_0)./(alpha*(e+e_2));% cycles/degree

res=pi*n*view_dist./180; % spatial resolution in cycles/degree

f_m=min(f_c,res./2);
hold on
h2=plot(f_m./view_dist./n.*180./pi)% cycles/degree to cycles per pixel

legend([h1,h2],{'Larcom','Wang'})


%% compare different fit functions

figure(6)
t=0:0.01:3;
a=exp(-t);
e2=2;
b=e2./(e2+3.*t);
c=1-t+t.*t./2-t.^3./6; % taylor series
% plot(t,[a;b;c]')
plot(t,[a;b]')



%% compare foveated images with 2 methods

data=FoveatedImageData(2);
image=data.image;

imageSize=size(image);

options1=struct();
options1.viewDist=1;
options1.numBands=20;
options1.filterHalfSize=20;
foveator1=foveateFilterBank(imageSize,data.fovea,options1);

options2=struct();
options2.numBands=20;
options2.filterHalfSize=20;
options2.cutoffMethod='Larcom';

options4=struct();
options4.numBands=20;
options4.filterHalfSize=20;
options4.cutoffMethod='Larcom';
n=min(imageSize);
options4.platSize=n/4;
options4.HPBW=n/4;
options4.minBandwidth = 0.1;

foveator2=foveateFilterBank(imageSize,data.fovea,options2);
foveator3=foveateSeperableOperator(imageSize,data.fovea, 20, 20);
foveator4=foveateFilterBank(imageSize,data.fovea,options4);
image1=foveator1.foveate(image);
image2=foveator2.foveate(image);
image3=foveator3.foveate(image);
image4=foveator4.foveate(image);

figure(7)
subplot(2,3,1)
imagesc(image1)
colormap gray
axis square off

subplot(2,3,2)
imagesc(image4)
% imagesc(image)
colormap gray
axis square off

subplot(2,3,3)
imagesc(image)
colormap gray
axis square off

subplot(2,3,4)
imagesc(image2)
colormap gray
axis square off

subplot(2,3,5)
imagesc(image3)
colormap gray
axis square off

subplot(2,3,6)
imagesc(image)
colormap gray
axis square off
