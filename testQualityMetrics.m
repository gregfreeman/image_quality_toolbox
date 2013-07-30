function results=testQualityMetrics(results,inputImage,outputImage)
% results=testQualityMetrics(results,inputImage,outputImage)
%
% results is a struct
% the output results struct will contain fields for different metrics of
% image quality
%
% inputImage - reference image
% outputImage - degraded image
%
% full reference:
% ms_ssim
% ssim
% mse
% relative_error
%
% no reference (blind):
% brisque
% bliinds2
% divine
%
% if results contains a 'settings' field with a struct containing a 'fovea'
% field, foveation based quality metrics will also be computed
% fssim
% fvmse + variants
%
%
%
%

%bound image pixels 0-255
boundedOutputImage=min(outputImage,255.*ones(size(outputImage)));
boundedOutputImage=max(boundedOutputImage,zeros(size(outputImage)));

results.ssim=ssim_index_m(inputImage,outputImage); 
[n,m]=size(inputImage);
if min(n,m) <= 128
    % adjust filter for ms_ssim to accomodate small images
    results.ms_ssim=msssim(inputImage,outputImage,[0.01 0.03],fspecial('gaussian', 7, 1.5));
    results.bounded_ms_ssim=msssim(inputImage,boundedOutputImage,[0.01 0.03],fspecial('gaussian', 7, 1.5));
else
    results.ms_ssim=msssim(inputImage,outputImage);
    results.bounded_ms_ssim=msssim(inputImage,boundedOutputImage);
end
results.relative_error=norm(inputImage(:)-outputImage(:))/norm(inputImage(:)); 
sqerror2=(inputImage-outputImage).^2;
results.mse=mean(sqerror2(:));

results.bounded_ssim=ssim_index_m(inputImage,boundedOutputImage); 
sqerror2=(inputImage-boundedOutputImage).^2;
results.bounded_mse=mean(sqerror2(:));
results.bounded_relative_error=norm(inputImage(:)-boundedOutputImage(:))/norm(inputImage(:)); 


if  isfield(results,'settings') && isfield(results.settings,'fovea')
%     fovea2=results.settings.fovea./size(inputImage); % fovea position relative to image frame
%     results.fvssim=fsbssim_index(inputImage,outputImage,fovea2); 
    results.fssim=fssim(inputImage,outputImage,struct('fovea',results.settings.fovea,'viewDist',3));
    if isfield(results.settings,'foveation_implementation')
        [results.fvmse,results.fvmse2,results.fvmse3]=image_mse_foveated(inputImage,outputImage,results.settings.fovea,results.settings.foveation_implementation); 
    else
        warning('foveation_implementation is not available');
    end
    %results.fvmse  % mse (f (reference), f(output) )
    %results.fvmse2 % mse (f (reference), output )
    %results.fvmse3 % mse (reference , output ) weighted by cutoff
    %frequency
end

% add dmos mappings for mse,ssim,ms_ssim so that they are on same scale
% as blind indices
try    
    fullpath = mfilename('fullpath');
    idx      = find(fullpath == filesep);
    thisfolder = fullpath(1:(idx(end)-1));
    dmos_fit=load(fullfile(thisfolder,'dmos_fit_fns'));
    % uses LIVE database with logistic regression
    results.mse_dmos=dmos_fit.fns.mse(results.mse); 
    % uses LIVE database with constrained polynomial regression
    results.ssim_dmos=dmos_fit.fns.ssim(results.ssim); % 
    % uses LIVE database with logistic regression
    results.ms_ssim_dmos=dmos_fit.fns.ms_ssim(results.ms_ssim);
catch exception
    results.mse_dmos=nan;
    results.ssim_dmos=nan;
    results.ms_ssim_dmos=nan;
    disp '*********** Error: '
    disp (getReport(exception,'extended'))
end
    
%brisque blind index
try
    results.brisque=brisque.brisquescore(outputImage);
catch exception
    results.brisque=nan;
    disp '*********** Error: '
    disp (getReport(exception,'extended'))
end
%bliinds2 blind index
try
    results.bliinds2=bliinds2.bliinds2_score(outputImage);
catch exception
    results.bliinds2=nan;
    disp '*********** Error: '
    disp (getReport(exception,'extended'))
end
%divine blind index
try
    results.divine=divine.divine(outputImage);
catch exception
    results.divine=nan;
    disp '*********** Error: '
    disp (getReport(exception,'extended'))
end

