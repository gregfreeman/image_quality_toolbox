 
% testQualityMetricsTest  test the testQualityMetrics function with a
% simple case
%
inputImage=double(imread('reference.png'));
outputImage=double(imread('degraded.png'));
results=struct();
results=testQualityMetrics(results,inputImage,outputImage);
disp(results)