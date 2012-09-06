% self=foveateSeperableOperator(imgSize,foveaLocations, numBands, filterHalfSize)
%    creates a fovation operator with precomputed filter banks and masks
%    for specific  fovea fixation locations.
%  imgSize - a 1x2 vector of the image size
%  foveaLocation - an nx2 vector of the location each of the fovea or
%       fixation points in [row column] form
%  numBands - the number of spatial filter banks to use when computing 
%       the foveated image
%  filterHalfSize - 1/2 the size (pixels) of the filter banks filters
%
%  returns an object with the foveate method.
%   call the foveate method with the actual image
%   i.e.  foveatedImage=self.fovate(image);
%
%  this operator uses the separability of the gaussian kernel for efficient
%  computation
%
% the foveation masks are computed from a half power roll off distance and
% a plateau size. 
% The two parameters are related to the viewing distance and image size

function self=foveateSeperableOperator(imgSize,foveaLocations, numBands, filterHalfSize)

if(nargin<2)
    foveaLocations=imgSize./2;
end
if(nargin<3)
    numBands=floor(log2(imgSize(1)))-5;
    if(numBands<2)
        numBands=2;
    end
end
if(nargin<4)
    filterHalfSize=20;
end
self.imgSize=imgSize;
self.foveaLocations=foveaLocations;
self.numBands=numBands;
self.filterHalfSize=filterHalfSize;

self.foveate=@foveate;
filterBanks=setupFilters(numBands, filterHalfSize);
masks=setupMasks(imgSize,foveaLocations);

    function masks=setupMasks(imgSize,foveaLocations)
        %size of the plateau in the foveation region as a pixel radius. 
        %platSize = 32;
        platSize = imgSize(2)/8;
        %minBandwidth = 0.2;
        minBandwidth = 0.1;
        %HPBW = 10;
        % half power bandwidth as a distance in pixels.  
        % The cutoff frequency rolls off exponentially as eccentricity increases. 
        % the HPBW distance is the number of pixels over which the cutoff
        % frequency halfs.
        HPBW = imgSize(1)/16;
        f_cutoff = repmat(minBandwidth, imgSize);  % cutoff frequency (relative scaled by Fs/2)
        % Fs is the pixel spatial freq
        lambda = log(2)/HPBW;
        rval = @(x,y,x0,y0) sqrt((x-x0).^2+(y-y0).^2)-platSize;

        %Place all of the foveae
        [c,r] = meshgrid(1:imgSize(2),1:imgSize(1));
        for iLocation=1:size(foveaLocations,1)
            %Note: x=col, y = row;
            % distance of each point to the fixation point (minus plateau radius)
            d = rval(c,r,foveaLocations(iLocation,2),foveaLocations(iLocation,1));
            d(d <0)=0;
            f_cutoff_i= exp(-d*lambda);

            %Use this number in areas where it exceeds the existing bandwidth
            f_cutoff= max(f_cutoff_i,f_cutoff);
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
        bandStep = 1/(numBands);
        fullSize = filterHalfSize*2 - 1;

        filterBanks = zeros([fullSize (numBands-1)]);
        cutoffs = NaN(numBands-1, 1);

        for iBand = 1:numBands
            cutoffs(iBand) = iBand * bandStep;
            sigma = sqrt(-2*log(0.5)) / (pi*cutoffs(iBand));
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
        
        for iBand=1:numBands-1
            f= conv2(filterBanks(:,iBand),filterBanks(:,iBand),padded,'valid');        
            foveated(masks(:,:,iBand))=f(masks(:,:,iBand)); 
            %figure(2), imagesc(foveated), colorbar
        end
        % no filtering for center of fovea
        foveated(masks(:,:,numBands))=image(masks(:,:,numBands)); 

    end
end



