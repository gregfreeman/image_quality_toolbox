% overall_mssim = fssim(img1, img2, options)
%    computes the foveated structural similarity index
% 
% 
% 
function [overall_fssim,fssim_map,debug_data] = fssim(img1, img2, options)

% Foveated Structural Similarity Index (F-SSIM)
%
%
%
% Z. Wang, E. P. Simoncelli and A. C. Bovik, "Multi-scale structural similarity
% for image quality assessment," Invited Paper, IEEE Asilomar Conference on
% Signals, Systems and Computers, Nov. 2003

if (nargin < 2 || nargin > 3)
   overall_ssim = -Inf;
   return;
end

if(~exist('options','var') || isempty(options))
    options=struct();
end
if(~isfield(options,'K'))
    options.K = [0.01 0.03]; 
end

if(~isfield(options,'win'))
   options.win = fspecial('gaussian', 11, 1.5);
end

if(~isfield(options,'levels'))
   options.levels = 5;
end

if(~isfield(options,'weight'))
   options.weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
end

if(~isfield(options,'method'))
   options.method = 'product';
end

[M N] = size(img1);
if(~isfield(options,'fovea'))
   options.fovea = [floor(M/2) floor(N/2)];
end

if(~isfield(options,'viewDist'))
   options.viewDist = 3;
end

K=options.K;
win=options.win;
levels=options.levels;
weight=options.weight;
method=options.method;
fovea=options.fovea;
viewDist=options.viewDist;

imgSize=size(img1);

if (size(img1) ~= size(img2))
   overall_mssim = -Inf;
   return;
end

[M N] = size(img1);
if ((M < 11) || (N < 11))
   overall_mssim = -Inf;
   return
end

if (length(K) ~= 2)
   overall_mssim = -Inf;
   return;
end

if (K(1) < 0 || K(2) < 0)
   overall_mssim = -Inf;
   return;
end
  
[H W] = size(win);

if ((H*W)<4 || (H>M) || (W>N))
   overall_mssim = -Inf;
   return;
end
   
if (levels < 1)
   overall_mssim = -Inf;
   return
end


min_img_width = min(M, N)/(2^(levels-1));
max_win_width = max(H, W);
if (min_img_width < max_win_width)
   overall_mssim = -Inf;
   return;
end

if (length(weight) ~= levels || sum(weight) == 0)
   overall_mssim = -Inf;
   return;
end

if ~((strcmp(method,'wtd_sum') || strcmp(method,'product')))
   overall_mssim = -Inf;
   return;
end

im1 = double(img1);
im2 = double(img2);

%levels = 5;
k_band = 5;
s = 0.5;
sep = 0.41;

bands1=subbands(im1,levels,sep,k_band,s);
bands2=subbands(im2,levels,sep,k_band,s);

normalize=0;
f=zeros(levels,1);
pixelWidth = 1/max(imgSize);
alpha = 0.106;
for iLevel = 1:levels
    f(iLevel) = s/(sep^-(iLevel-1))*0.0175*viewDist/pixelWidth;
    normalize = max(normalize,exp(-alpha*f(iLevel)));
end

fssim_map=ones(imgSize-[10 10]);
fssim_band=zeros([(imgSize-[10 10]) levels]);
ssim_array=zeros(levels,1);
ssim_map_array=cell(levels,1);
cs_array=zeros(levels,1);
cs_map_array=cell(levels,1);
Sf_array=zeros([(imgSize-[10 10]) levels]);
debug_data=cell(3,1);
debug_data{1}=cell(levels,1);
debug_data{2}=cell(levels,1);
debug_data{3}=cell(levels,1);

for iLevel = 1:levels
   [ssim_array(iLevel) ssim_map_array{iLevel} cs_array(iLevel) cs_map_array{iLevel}] = ...
       ssim_index_new(bands1{iLevel}, bands2{iLevel}, K, win);
   debug_data{1}{iLevel}= ssim_map_array{iLevel};
   [Sf,e]=foveateSensitivity(imgSize,fovea, f(iLevel), viewDist);
   Sf=Sf./normalize;
%    Sf(f(iLevel)>f_m)=0;
   mapSize=size(ssim_map_array{iLevel});
   offset=(imgSize-mapSize)./2;
   SfTrim=Sf(offset(1)+1:offset(1)+mapSize(1), ...
       offset(2)+1:offset(2)+mapSize(2));
   debug_data{2}{iLevel}= SfTrim;
   map=ssim_map_array{iLevel};
   map(map<0)=0;
   Sf_array(:,:,iLevel)=SfTrim;
   fssim_band(:,:,iLevel)=map.^SfTrim;
   debug_data{3}{iLevel}= fssim_band(:,:,iLevel);
   fssim_map=fssim_map.*fssim_band(:,:,iLevel);
   debug_data{4}=e;
end

overall_fssim =mean2(fssim_map);


function bands=subbands(img,levels,sep,k,s)
bands=cell(levels,1);
sp = sqrt(log(k))/(pi*s*sqrt(k^2-1));
% hsize=[11 11];
for iLevel=1:(levels-1)
    sigma1=sp*sep^(-(iLevel-1));
    sigma2=sp*k*sep^(-(iLevel-1));
%     hsize=ceil(sigma1*3);
    ksize=floor(sigma2*8+1+0.5);
%     if(ksize>hsize(1))
%         warning('small kernel size')
%     end

    kernel1=fspecial('gaussian', [ksize ksize],sigma1);
    kernel2=fspecial('gaussian', [ksize ksize],sigma2);

    m1=imfilter(img,kernel1,'symmetric','same');
    m2=imfilter(img,kernel2,'symmetric','same');
    bands{iLevel}=m2-m1;
end
sigma1=s*sep^(-(levels-1));
ksize=floor(sigma1*8+1+0.5);
% if(ksize>hsize(1))
%     warning('small kernel size')
% end
kernel1=fspecial('gaussian', [ksize ksize],sigma1);
bands{levels}=imfilter(img,kernel1);
