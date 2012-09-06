function [mse,mse2,fmse]=image_mse_foveated(inputImage,outputImage,fovea,foveator)
%[mse,mse2,fmse]=image_mse_foveated(inputImage,outputImage,fovea,foveator)
%  mse is the mse between the foveated reference and foveated output
%  mse2 is the mse between the foveated reference and unaltered output
%  fmse is the mse between the unaltered reference and unaltered output 
%    weighted by the foveation cutoff frequency
%

[n1,n2]=size(inputImage);
if ~exist('foveator','var')
    options1=struct();
    options1.viewDist=1;
    options1.numBands=20;
    options1.filterHalfSize=20;
    foveator=foveateFilterBank([n1,n2],fovea,options1);
end
foveatedInputImage=foveator.foveate(inputImage);
foveatedOutputImage=foveator.foveate(outputImage);

sqerror=(foveatedInputImage-foveatedOutputImage).^2;
mse=mean(sqerror(:));
sqerror2=(foveatedInputImage-outputImage).^2;
mse2=mean(sqerror2(:));

f_cutoff=foveateCutoffFreq([n1,n2],fovea, foveator.options);
f_cutoff_sum=sum(f_cutoff(:));
sqerror3=(inputImage-outputImage).^2;
sqerror3_weighted=sqerror3.*f_cutoff.^2./f_cutoff_sum;
fmse=sum(sqerror3_weighted(:));
