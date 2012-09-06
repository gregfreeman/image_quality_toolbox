%testing for fsbssim_index 

addpath('C:/Development/Contracts/RFCode/branches/diagnostics/cpp/SB-SSIM/Release')
addpath('C:/Development/Contracts/RFCode/branches/diagnostics/cpp/SB-SSIM/Debug')
addpath('C:/Development/Contracts/RFCode/branches/diagnostics/cpp/SB-SSIM/MexTest')

fov = [180/512 256/512];

nTests=2;
q=zeros(1,nTests);
Qarray=cell(1,nTests);
t=zeros(1,nTests);
k=5;
sep=0.41;
l=5;
weights=[0.0448,0.2856,0.3001,0.2363,0.1333];

fName ='Cameraman.BMP';
image1= imread(fName);
fName2 = 'Cameraman_MS0.BMP';
image2=imread(fName2);

% offset=128+64;
% image1s=image1(offset+1:offset+128,offset+1:offset+128);
% image2s=image2(offset+1:offset+128,offset+1:offset+128);

iCase = 1;
%refs={ image1, double(image1) , uint16(image1)  } ;
%tests={ image2, double(image2), image1 } ;
for refs={ double(image1) }
%for iRef=1:2
%    for iTest=1:2
%     for tests={ double(image2s), double(image1s) }
    for tests={ double(image2) }
    %    ref=refs{iRef};
    %    test=tests{iTest};
        test=tests{1};
        ref=refs{1};

        whos test
        whos ref
        
%         try
            display fsbssim_index
            tic
            imgSize=size(test);
            [q(iCase) ,Qarray{iCase},debug_data_mex]=fsbssim_indexd(test,ref,fov);
            t(iCase)=toc;
            iCase=iCase+1;
%         catch exception
%             disp(exception)
%             disp(exception.stack)
%         end
%                
                
%         try
            display fssim
            tic
            [q(iCase) ,Qarray{iCase},debug_data_m]=fssim(test,ref,struct('fovea',fov.*imgSize));
            t(iCase)=toc;
            iCase=iCase+1;
%         catch exception
%             disp(exception)
%         end
               
    end
end

imagesc([Qarray{1} Qarray{2}.*255])
colormap gray
colorbar


debug_data_mex=debug_data_mex([3,1,2,4]);

figure(1)
for i=1:5
    subplot(3,1,1)
    imagesc([debug_data_mex{1}{i} debug_data_m{1}{i}])
    colormap gray
    colorbar
    subplot(3,1,2)
    imagesc([debug_data_mex{2}{i} debug_data_m{2}{i}])
    colormap gray
    colorbar
    subplot(3,1,3)
    imagesc([debug_data_mex{3}{i} debug_data_m{3}{i}])
    colormap gray
    colorbar
    %pause(0.5)
end

figure(2)
e=debug_data_m{4}(6:507,6:507);
imagesc([debug_data_mex{4} e])
colormap gray
colorbar

% stack image levels together
sz=502;
a=zeros(sz*5,sz*2);
b=zeros(sz*5,sz*2);
c=zeros(sz*5,sz*2);
for i=1:5
    a(1+(i-1)*sz:i*sz,:)=[debug_data_mex{1}{i} debug_data_m{1}{i}];
    b(1+(i-1)*sz:i*sz,:)=[debug_data_mex{2}{i} debug_data_m{2}{i}];
    c(1+(i-1)*sz:i*sz,:)=[debug_data_mex{3}{i} debug_data_m{3}{i}];
end
figure(3)
imagesc(a)
colormap gray
colorbar
axis off
axis square
% figure(4)
% imagesc(b)
% colormap gray
% colorbar
% axis off
% axis square
figure(5)
imagesc(c)
colormap gray
colorbar
axis off
axis square
%pause(0.5)

figure(6)
imagesc([Qarray{1} Qarray{2}.*255])
colormap gray
colorbar
axis off
axis square
