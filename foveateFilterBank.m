% self=foveateFilterBank(imgSize,foveaLocations, options)
%    creates a fovation operator with precomputed filter banks and masks
%    for specific  fovea fixation locations.
%  imgSize - a 1x2 vector of the image size
%  foveaLocation - an nx2 vector of the location each of the fovea or
%       fixation points in [row column] form
%  options - struct for configuring options
%  options.numBands - the number of spatial filter banks to use when computing 
%       the foveated image.  default is a function of the size of the image
%  options.filterHalfSize - 1/2 the size (pixels) of the filter banks filters
%     defaults to 20
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
%  returns an object with the foveate method.
%   call the foveate method with the actual image
%   i.e.  foveatedImage=self.fovate(image);
%
%  this operator uses the separability of the gaussian kernel for efficient
%  computation, but computes each filter across the entire image. (so could
%  be made more efficient)
%

function self=foveateFilterBank(imgSize,foveaLocations, options)

if(nargin<2)
    foveaLocations=imgSize./2;
end
if(nargin<3)
    options=struct();
end

if ~isfield(options,'numBands')
    options.numBands =floor(log2(imgSize(1)))-5;
    if(options.numBands <2)
        options.numBands =2;
    end
end
numBands=options.numBands;

if ~isfield(options,'filterHalfSize')
    options.filterHalfSize=20;
end
filterHalfSize=options.filterHalfSize;

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
self.filterBanks=setupFilters(options.numBands, options.filterHalfSize);
self.masks=setupMasks(imgSize,foveaLocations,options);
self.options=options;

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

    function masks=setupMasks(imgSize,foveaLocations,options)
        % create a bank masks of which correspond to which filter to use
        % for different locations within the image
        
        if strcmp(options.cutoffMethod,'Larcom')
            f_cutoff=cutoffFreqFromHPBWAndPlateau(imgSize,foveaLocations,options.HPBW,options.platSize,options.minBandwidth);
        else
            f_cutoff=cutoffFreqFromViewDist(imgSize,foveaLocations,options.viewDist);
        end
                
        % create binary masks for which filter corresponds with each point in
        % the image
        %
        % mask index is an index of the filter bank for each point 
        % band 
        maskIndex=ceil(f_cutoff*numBands);
        maskIndex(maskIndex<1)=1;
        maskIndex(maskIndex>numBands)=numBands;
        
        masks=false([imgSize numBands]);
        for iBand=1:numBands
            masks(:,:,iBand)=maskIndex==iBand;
        end
    end

    function filterBanks=setupFilters(numBands, filterHalfSize)
        % create a bank of low pass filters
        bandStep = 1/(numBands);
        fullSize = filterHalfSize*2 - 1;

        filterBanks = zeros([fullSize (numBands-1)]);
        cutoffs = NaN(numBands, 1);

        for iBand = 1:numBands
            cutoffs(iBand) = iBand * bandStep;
            sigma = sqrt(2*log(2)) / (pi*cutoffs(iBand));
            f=fspecial('gaussian', fullSize, sigma);
            stencil = rot90(f,2);
            [u,s,v] = svd(stencil);
            hcol = u(:,1) * sqrt(s(1));
            hrow = conj(v(:,1)) * sqrt(s(1));
            if(hcol~=hrow)
                disp('Why are the vectors not the same')
            end

            filterBanks(:,iBand) = hcol;
        end
        
    end

    function foveated=foveate(image)
        
        foveated=zeros(size(image));
        % pad image so that after applying filters the output will be the
        % same size
        filterLength=filterHalfSize-1;
        padded=zeros(size(image) + [1 1].*filterLength.*2);
        padded(filterLength+1:end-filterLength,filterLength+1:end-filterLength)=image;
        % corners
        padded(1:filterLength        ,1:filterLength        )=image(1  ,1  );
        padded(1:filterLength        ,end-filterLength+1:end)=image(1  ,end);
        padded(end-filterLength+1:end,1:filterLength        )=image(end,1  );
        padded(end-filterLength+1:end,end-filterLength+1:end)=image(end,end);
        
        % sides
            %top/bottom
        padded(1:filterLength        ,filterLength+1:end-filterLength )=repmat(image(1,  :),[filterLength 1]);
        padded(end-filterLength+1:end,filterLength+1:end-filterLength )=repmat(image(end,:),[filterLength 1]);
            %left/right
        padded(filterLength+1:end-filterLength ,1:filterLength        )=repmat(image(:,1  ),[1 filterLength]);
        padded(filterLength+1:end-filterLength ,end-filterLength+1:end)=repmat(image(:,end),[1 filterLength]);      

        % filter image with a bank of filters
        % mask the filter output for the appropriate image location
        for iBand=1:numBands-1
            f= conv2(self.filterBanks(:,iBand),self.filterBanks(:,iBand),padded,'valid');        
            foveated(self.masks(:,:,iBand))=f(self.masks(:,:,iBand)); 
        end
        % no filtering for center of fovea
        foveated(self.masks(:,:,numBands))=image(self.masks(:,:,numBands)); 

    end
end



