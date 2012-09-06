% w=foveateSensitivity(imgSize,foveaLocations, f, view_dist)
%    creates a fovation weigh matrix for an image
%        given specific fovea fixation locations.
%  imgSize - a 1x2 vector of the image size 
%  foveaLocation - an nx2 vector of the location each of the fovea or
%       fixation points in [row column] form
%  f - spatial frequency in cycles/degree
%  view_dist - viewing distance as a mulitple of screen size
%
% [1] Z. Wang and a C. Bovik, “Embedded foveation image coding,” IEEE
% Transactions on Image Processing, vol. 10, no. 10, pp. 1397-410.
%
function [Sf,e]=foveateSensitivity(imgSize,foveaLocations, f, view_dist)

n1=imgSize(1);
n2=imgSize(2);

if(~exist('foveaLocations','var') || isempty(foveaLocations))
    foveaLocations=imgSize./2;
end

if(~exist('f','var') || isempty(f))
    fmax=pi*n1*view_dist./360; % spatial resolution in cycles/degree
    fmin=pi*1*view_dist./360; % spatial resolution in cycles/degree    
    f=sqrt(fmax*fmin);% geometric mean
end

if(~exist('view_dist','var')  || isempty(view_dist))
    view_dist=1;
end

alpha=0.106;
e_2=2.3;
CT_0=1/64;

distance = @(x,y,x0,y0) sqrt((x-x0).^2+(y-y0).^2);
d=Inf*ones(n1,n2); % distance from fixation points
% indices for image
[x2,x1] = meshgrid(1:n1,1:n2);
for iLocation=1:size(foveaLocations,1)
    %Note: x2=col, x1 = row;
    d = min(d,distance(x1,x2,foveaLocations(iLocation,1),foveaLocations(iLocation,2)));
end

e=atan(d./view_dist./n1).*180./pi; % eccentricity in degrees.
res=pi*n1*view_dist./180; % spatial resolution in cycles/degree

% f[i] = s/pow(sep,-i)*0.0175*viewDistance/pixelWidth;
% norm = max(norm,exp(-alpha*f[i]));
% e[y][x] = min(e[y][x],atan(sqrt(pow(fovea[i]*a.cols-x,2)+pow(fovea[numFovea+i]*a.rows-y,2))*pixelWidth/viewDistance)*180.0/PI);			
% double w = exp(-alpha*f[i]*(e[y+5][x+5]+e2)/e2)/norm;
% e_prime = atan(d*pixelWidth/viewDistance)*180.0/pi;			

Sf_prime = exp(-alpha.*f.*(e+e_2)./e_2);

%CT=CT_0*exp(alpha.*f.*(e+e_2)./e_2);
%CS=1./CT;
f_c=e_2*log(1/CT_0)./(alpha*(e+e_2));% cycles/degree
% f_d=res./2;
% f_m=min(f_c,f_d);
Sf=exp(-0.0461*f.*(e+e_2));
Sf(f>f_c)=0;

