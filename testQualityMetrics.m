function results=testQualityMetrics(results,inputImage,outputImage)

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
results.relative_error=norm(inputImage-outputImage)/norm(inputImage); 
sqerror2=(inputImage-outputImage).^2;
results.mse=mean(sqerror2(:));

results.bounded_ssim=ssim_index_m(inputImage,boundedOutputImage); 
sqerror2=(inputImage-boundedOutputImage).^2;
results.bounded_mse=mean(sqerror2(:));
results.bounded_relative_error=norm(inputImage-boundedOutputImage)/norm(inputImage); 


if isfield(results.settings,'fovea')
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
