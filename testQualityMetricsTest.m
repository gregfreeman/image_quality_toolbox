


inputImage=double(imread('reference.png'));
outputImage=double(imread('degraded.png'));
results=struct();
results=testQualityMetrics(results,inputImage,outputImage);
