function qualityscore  = brisquescore(imdist)

import brisque.*

fullpath = mfilename('fullpath');
idx      = find(fullpath == filesep);
bfolder = fullpath(1:(idx(end)-1));

if(size(imdist,3)==3)
    imdist = uint8(imdist);
    imdist = rgb2gray(imdist);
end

imdist = double(imdist);

if(nargin<2)
feat = brisque_feature(imdist);
% disp('feat computed')
end


%---------------------------------------------------------------------
%Quality Score Computation
%---------------------------------------------------------------------

if 0
    fid = fopen(fullfile(bfolder,'test_ind'),'w');

    for jj = 1:size(feat,1)

    fprintf(fid,'1 ');
    for kk = 1:size(feat,2)
    fprintf(fid,'%d:%f ',kk,feat(jj,kk));
    end
    fprintf(fid,'\n');
    end

    fclose(fid);

    % warning off all
    % delete output test_ind_scaled dump

    delete (fullfile(bfolder,'output'))
    delete (fullfile(bfolder,'test_ind_scaled'))
    delete (fullfile(bfolder,'dump'))

    % system('svm-scale -r allrange test_ind >> test_ind_scaled');
    cmd=sprintf('"%s/svm-scale" -r "%s/allrange" "%s/test_ind" >> "%s/test_ind_scaled"',bfolder,bfolder,bfolder,bfolder);
    disp(cmd)
    system(cmd)

    % system('svm-predict -b 1 test_ind_scaled allmodel output >>dump');
    cmd=sprintf('"%s/svm-predict" -b 1 "%s/test_ind_scaled" "%s/allmodel" "%s/output"',bfolder,bfolder,bfolder,bfolder);
    disp(cmd)
    system(cmd)

    load(fullfile(bfolder,'output'))
else
    % use clone git://github.com/gregfreeman/libsvm.git branch -b new_matlab_interface
    [range,lu]=libsvmreadrange(fullfile(bfolder,'allrange'));
    sfeat=libsvmscale(feat,range,lu);
    allmodel=libsvmreadmodel(fullfile(bfolder,'allmodel'),36);
    [output, p ] = svmpredict(sfeat, allmodel);
%     disp(output)
%     disp(p)
end

qualityscore = output;
