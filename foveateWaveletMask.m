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

function self=foveateWaveletMask(imgSize,foveaLocations, numBands, filterHalfSize)

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

self.getMask=@getMask;
filterBanks=setupFilters(numBands, filterHalfSize);
masks=setupMasks(imgSize,foveaLocations);

    function masks=setupMasks(imgSize,foveaLocations)
        %platSize = 32;
        platSize = imgSize(2)/8;
        %minBandwidth = 0.2;
        minBandwidth = 0.1;
        %HPBW = 10;
        HPBW = imgSize(1)/16;
        maskBins = repmat(minBandwidth, imgSize);
        lambda = log(2)/HPBW;
        rval = @(x,y,x0,y0) sqrt((x-x0).^2+(y-y0).^2)-platSize;

        %Place all of the foveae
        [c,r] = meshgrid(1:imgSize(2),1:imgSize(1));
        for iLocation=1:size(foveaLocations,1)
            %Note: x=col, y = row;
            curr = rval(c,r,foveaLocations(iLocation,2),foveaLocations(iLocation,1));
            fovControl = exp(-curr*lambda);
            fovControl(curr<0) = 1;

            %Use this number in areas where it exceeds the existing bandwidth
            maskBins(fovControl>maskBins) = fovControl(fovControl>maskBins);
        end
        maskIndex=ceil(maskBins*numBands);
        maskIndex(maskIndex<1)=1;
        maskIndex(maskIndex>numBands)=numBands;
        masks=false([imgSize numBands]);
        %figure(1), imagesc(maskIndex), colorbar
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
   

    function mask=getMask
        
        mask=zeros(self.imgSize);
        mask(1,1)=1;
        
        for k=1:(log2(self.imgSize(1))-1)
            k2=2.^k;
        end
        % no filtering for center of fovea
%         foveated(masks(:,:,numBands))=image(masks(:,:,numBands)); 

    end
end



