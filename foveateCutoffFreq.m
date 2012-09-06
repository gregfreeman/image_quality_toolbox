% f_cutoff=foveateCutoffFreq(imgSize,foveaLocations, options)
%   Finds the cutoff frequency across a foveated image
%   returns array of cutoff frequencies
%       frequencies are relative to Fs/2 or pi radians
%
%   Finds the cutoff frequency across a foveated image
%  imgSize - a 1x2 vector of the image size
%  foveaLocation - an nx2 vector of the location each of the fovea or
%       fixation points in [row column] form
%  options - struct for configuring options
%  options.cutoffMethod - function for computing the local cutoff frequency
%  across the image.  Parameter is a string
%   'Wang'  (default)
%     Z. Wang and A. C. Bovik, “Embedded foveation image coding,” IEEE
%     Transactions on Image Processing, vol. 10, no. 10, pp. 1397–410, Oct.
%     2001.
%   'Larcom'
%     R. Larcom and T. R. Coffman, “Foveated image formation through
%     compressive sensing,” 2010 IEEE Southwest Symposium on Image Analysis
%     & Interpretation (SSIAI), pp. 145–148, 2010.
%
%  'Wang' options:
%   options.viewDist (defaults to 3)
%         ratio of the viewer distance to the screen size
%
%   'Larcom' options:
%    options.HPBW   (defaults to imgSize(1)/16)
%         half power bandwidth as a distance in pixels.  
%         The cutoff frequency rolls off exponentially as eccentricity increases. 
%         the HPBW distance is the number of pixels over which the cutoff
%         frequency halfs.
%    options.platSize (defaults to = imgSize(2)/8)
%         size of the plateau in the foveation region as a pixel radius.
%    options.minBandwidth   (defaults to 0.1)
%         minBandwidth is a minimum on the cutoff frequency (relative to
%         Fs/2)
%

function f_cutoff=foveateCutoffFreq(imgSize,foveaLocations, options)

if(nargin<2)
    foveaLocations=imgSize./2;
end
if(nargin<3)
    options=struct();
end


self.imgSize=imgSize;
self.foveaLocations=foveaLocations;

self.foveate=@foveate;

if ~isfield(options,'cutoffMethod')
    options.cutoffMethod='Wang';
end

if strcmp(options.cutoffMethod,'Wang')
%     Z. Wang and A. C. Bovik, “Embedded foveation image coding,” IEEE
%     Transactions on Image Processing, vol. 10, no. 10, pp. 1397–410, Oct.
%     2001.

% view distance is a ratio of the distance from a screen to the viewer 
    %  and the image size
    if ~isfield(options,'viewDist')
        options.viewDist=3;
    end

elseif strcmp(options.cutoffMethod,'Larcom')
    
%     R. Larcom and T. R. Coffman, “Foveated image formation through
%     compressive sensing,” 2010 IEEE Southwest Symposium on Image Analysis
%     & Interpretation (SSIAI), pp. 145–148, 2010.

    % half power bandwidth as a distance in pixels.  
    % The cutoff frequency rolls off exponentially as eccentricity increases. 
    % the HPBW distance is the number of pixels over which the cutoff
    % frequency halfs.
    if ~isfield(options,'HPBW')
        options.HPBW= imgSize(1)/16;
    end
    
    %size of the plateau in the foveation region as a pixel radius.
    if ~isfield(options,'platSize')
        options.platSize= imgSize(2)/8;
    end
    if ~isfield(options,'minBandwidth')
        options.minBandwidth=0.1;
    end

else
    warning('Unrecognized foveation method')
end


    % create a bank masks of which correspond to which filter to use
    % for different locations within the image

    if strcmp(options.cutoffMethod,'Larcom')
        f_cutoff=cutoffFreqFromHPBWAndPlateau(imgSize,foveaLocations,options.HPBW,options.platSize,options.minBandwidth);
    else
        f_cutoff=cutoffFreqFromViewDist(imgSize,foveaLocations,options.viewDist);
    end
                


    % cutoffFreqFromViewDist(imgSize,foveaLocations,viewDist)
    % implementation matches model from Wang's work
    % returns frequency (scaled by Fs/2)
    % view distance is a ratio of the distance from a screen to the viewer 
    %  and the image size
    function f_cutoff=cutoffFreqFromViewDist(imgSize,foveaLocations,viewDist)
        alpha=0.106;
        e_2=2.3;
        CT_0=1/64;

        distance = @(x,y,x0,y0) sqrt((x-x0).^2+(y-y0).^2);
        d=Inf*ones(imgSize); % distance from fixation points

        %Place all of the foveae
        [c,r] = meshgrid(1:imgSize(2),1:imgSize(1));
        for iLocation=1:size(foveaLocations,1)
            %Note: x=col, y = row;
            % distance of each point to the fixation point (minus plateau radius)
            d1=distance(c,r,foveaLocations(iLocation,2),foveaLocations(iLocation,1));
            d = min(d,d1);
        end
        
        n1=min(imgSize);
        e=atan(d./viewDist./n1).*180./pi; % eccentricity in degrees.
        res=pi*n1*viewDist./180; % spatial resolution in cycles/degree
        f_c=e_2*log(1/CT_0)./(alpha*(e+e_2));% cycles/degree

        f_cutoff= min(f_c./res,1);% scaled frequency
    end

    %cutoffFreqFromHPBWAndPlateau(imgSize,foveaLocations,HPBW,platSize,minBandwidth)
    % implementation matches model from Ron Larcom's work
    % returns frequency (scaled by Fs/2)
    % Fs is the pixel spatial freq
    % half power bandwidth as a distance in pixels.  
    % The cutoff frequency rolls off exponentially as eccentricity increases. 
    % the HPBW distance is the number of pixels over which the cutoff
    % frequency halfs.
    % platSize=size of the plateau in the foveation region as a pixel radius. 
    % minBandwidth is a minimum on the cutoff frequency
    function f_cutoff=cutoffFreqFromHPBWAndPlateau(imgSize,foveaLocations,HPBW,platSize,minBandwidth)
        lambda = log(2)/HPBW;
        distance = @(x,y,x0,y0) sqrt((x-x0).^2+(y-y0).^2);
        d=Inf*ones(imgSize); % distance from fixation points

        %Place all of the foveae
        [c,r] = meshgrid(1:imgSize(2),1:imgSize(1));
        for iLocation=1:size(foveaLocations,1)
            %Note: x=col, y = row;
            % distance of each point to the fixation point (minus plateau radius)
            d1=distance(c,r,foveaLocations(iLocation,2),foveaLocations(iLocation,1))-platSize;
            d1(d1 <0)=0;
            d = min(d,d1);
        end
        f_cutoff= exp(-d*lambda);
        % set minimum frequency
        f_cutoff = max(f_cutoff,repmat(minBandwidth, imgSize));  % cutoff frequency (relative scaled by Fs/2)
    end

end



