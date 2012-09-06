% w=foveateWaveletSensitivity(imgSize,foveaLocations, levels, view_dist,
% options)
%    creates a fovation weigh matrix for wavelet coefficients 
%        given specific fovea fixation locations.
%  imgSize - a 1x2 vector of the image size (same as the wavelet
%  coefficient matrix
%  foveaLocation - an nx2 vector of the location each of the fovea or
%       fixation points in [row column] form
%  levels - the number of wavelet decomposition levels
%  view_dist - viewing distance as a mulitple of screen size
%  options - structure for variation options 
%    options.use_band_weighting - 0 or 1 on whether to include the subband
%    weighting term based on empirical noise experiments (defaults to 0)
%    options.plot_band_weighting - 0 or 1 on whether to plot the affects of
%     using the band weighting option
%
% [1] Z. Wang and a C. Bovik, “Embedded foveation image coding,” IEEE
% Transactions on Image Processing, vol. 10, no. 10, pp. 1397-410.
%
function w=foveateWaveletSensitivity(imgSize,foveaLocations, levels, view_dist, options)

if(nargin<2 || isempty(foveaLocations))
    foveaLocations=imgSize./2;
end
if(nargin<3 || isempty(levels))
    levels=floor(log2(imgSize(1)))-4;
    if(levels<2)
        levels=2;
    end
end
if(nargin<4 || isempty(view_dist))
    view_dist=1;
end
if(nargin<5 || isempty(options))
    options=struct();
end
if(~isfield(options,'use_band_weighting'))
    options.use_band_weighting=0; 
end
if(~isfield(options,'plot_band_weighting'))
    options.plot_band_weighting=0; 
end

alpha=0.106;
e_2=2.3;
CT_0=1/64;

%size of level 1 panel
n1=imgSize(1)/2;
n2=imgSize(2)/2;
distance = @(x,y,x0,y0) sqrt((x-x0).^2+(y-y0).^2);
d_level=Inf*ones(n1,n2); % distance from fixation points
% indices for a subset of wavelet coeficients
[x2,x1] = meshgrid(1:n1,1:n2);
for iLocation=1:size(foveaLocations,1)
    %Note: x2=col, x1 = row;
    d_level = min(d_level,distance(x1,x2,foveaLocations(iLocation,1)./2,foveaLocations(iLocation,2)./2));
end

% ---------
% | a | v |
% |-------|
% | h | d |
% ---------
%         h   v   d
aOffset=[1 0;0 1;1 1];

e_level=atan(d_level.*2./view_dist./n1)*180./pi; % eccentricity in degrees.
res=pi*n1*view_dist./180; % spatial resolution in cycles/degree

e=zeros(n1*2,n2*2);
f=zeros(n1*2,n2*2);
% fill eccentricity and spatial freqency of each point in wavelet domain
s1=n1;
s2=n2;
for iLevel=1:levels
    for iOrientation=1:3
        e(aOffset(iOrientation,1)*s1+1:(aOffset(iOrientation,1)+1)*s1, ...
          aOffset(iOrientation,2)*s2+1:(aOffset(iOrientation,2)+1)*s2) = e_level;
        f(aOffset(iOrientation,1)*s1+1:(aOffset(iOrientation,1)+1)*s1, ...
          aOffset(iOrientation,2)*s2+1:(aOffset(iOrientation,2)+1)*s2) = res*2^-iLevel;
    end
    if iLevel<levels % half size for next level (unless last iteration)
        s1=s1/2;
        s2=s2/2;
        e_level=e_level(1:2:end,1:2:end);    
    end
end
% LL_levels (scaling coefficients)
e(1:s1,1:s2)= e_level; 
f(1:s1,1:s2)= res*2^-levels;

%CT=CT_0*exp(alpha.*f.*(e+e_2)./e_2);
%CS=1./CT;
f_c=e_2*log(1/CT_0)./(alpha*(e+e_2));% cycles/degree
f_d=res./2;
f_m=min(f_c,f_d);
Sf=exp(-0.0461*f.*(e+e_2));
Sf(f>f_m)=0;


if ~options.use_band_weighting
    w=Sf;
else
    
    % sub-band weighting - details in [1]
    Beta1=1;
    Beta2=2.5;
%     a=0.495;
%     k=0.466;
%     f_0=0.401;
%     g_theta=[ 1.501 1 1 0.534];

    Sw=zeros(n1*2,n2*2);

    s1=n1;
    s2=n2;
    Sw_table=[.3842 .27 .1316; .3818 .3326 .2138; .2931 .3019 .2442; .1804 .2129 .2098; .0905 .1207 .1430; .0372 .0558 .0791];
    Sw_table=Sw_table(:,[1 2 2 3]);
    for iLevel=1:levels
        for iOrientation=1:3
            Sw_level=Sw_table(iLevel,iOrientation+1);
            Sw(aOffset(iOrientation,1)*s1+1:(aOffset(iOrientation,1)+1)*s1, ...
              aOffset(iOrientation,2)*s2+1:(aOffset(iOrientation,2)+1)*s2) = Sw_level;
        end
        if iLevel<levels
            s1=s1/2;
            s2=s2/2;
        end
    end
    Sw_level=Sw_table(levels,1);
    Sw(1:s1,1:s2)= Sw_level;
%    w=(Sw.^Beta1.*Sf.^Beta2); % as described in [1]
    w=(Sw.^(Beta1/Beta2).*Sf); % has closer scaling to alternative
end

% diagnostic to show affect of band_weighting
if options.plot_band_weighting && options.use_band_weighting
  figure(1)
  subplot(1,3,1)
  imagesc(Sf)
  colormap gray
  colorbar
  axis off
  axis square
  title 'Sf'
  subplot(1,3,2)
  imagesc(Sw)
  colormap gray
  colorbar
  axis off
  axis square
  title 'Sw'
  subplot(1,3,3)
  imagesc(w)
  colormap gray
  colorbar
  axis off
  axis square
  title 'w'
end




