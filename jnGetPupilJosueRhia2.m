%adjusted from fileexchange scripts
%JN 8/20/2018, modeled after Andrew's original getPupil script
%load bin file
% pathy='Z:\pupilMuscimol';

clear all
clear all

close all
close all

showIm = false;
saveVideo = false; 
savePupil=true;
%pupilVidName=[path '\m3_Pupil_T174846'];
pupilVidName=date;
date='m585L_Pupil_T133932'
filename = [date,'.bin'];

%%%%%%%%%
conthres=0.4;
%%%%%%%%%
%m585L need [10 20] and conthres 0.4 and small ROI
%mouse10 needs 0.4 contrast
%TIP: PICK a small ROI within the eye so that its easier to pick up pupil

fileInfo = dir(filename);
% dbin = dir(fullfile(currentFolder,inDir,'*.bin'));
% for ibin = 1:length(dbin)
%     if(contains(dbin(ibin).name,'pupil'))
%         pupilFile = fullfile(currentFolder,inDir,dbin(ibin).name);
%     end
% end

%adjusting:
fileInfo = dir(filename);
headerBytes = 24;
fileSize = fileInfo.bytes - headerBytes;
fid = fopen(filename);

header = fread(fid,headerBytes/8,'double'); % header is written as doubles (= 8 bytes)
w = header(1);
h = header(2);
b = header(3); % Should be using grayscale so doesn't matter. divide by 8 bytes

imSize = w*h;
nFrames = fileSize/(imSize+8);
%above adjusted

% below old style
% fileInfo = dir(pupilFile);
% fileSize = fileInfo.bytes;
% 
% fid = fopen(pupilFile);
% 
% header = fread(fid,3,'double');
% w = header(1);
% h = header(2);
% b = header(3)/8;
% imSize = w*h*b;
% 
% nFrames = (fileSize - 3*8) / (imSize + 8);
%above old style
if(mod(nFrames,1)~=0)
    error('Total size is not a multiple of image size... something is wrong');
end

if(saveVideo)
    F = getframe; %begin storing frames
end
pupil_diam = zeros(nFrames,1);
pupil_tstamp = zeros(nFrames,1);

prevDiam=15; %these checks were added because of issues in individual experiments (images to bright or too dark or eyelashes etc..)
%It would be best to get these set such that they don't need to be adjusted
%for every video. This depends partly on getting a system for avoiding
%light issues during recording. The mouse hat with foam works well when it
%works..but some mice can see the foam and pull off their hat. Perhaps make
%hats more stylish.
diamCheck=10; %max pixel change in diameter allowed between successive frames


%%%set the roi for the pupil video
    pupil_tstamp(1) = fread(fid,1,'double');
    %     d=fread(fid,imSize,'ubit8');
    %     dd=uint8(fread(fid,imSize,'ubit8'));
    %     dd=uint8(fread(fid,imSize,'ubit8'));
    %     data=dd;
    %data = uint8(fread(fid,imSize,'ubit8'));
    frame = uint8(fread(fid,imSize,'uint8'));
    Data = reshape(frame, [h,w]);
    Data = flip(transpose(Data));   
%     Data = flip(transpose(Data));  
%     Data = flip(transpose(Data));   
%     Data = imadjust(Data);


    [x2,y2,BW,xi2,yi2]=roipoly(Data);
   

endloopwith=nFrames-1;
% endloopwith=500;
ctr=0;
for iF = 1:endloopwith
    pupil_tstamp(iF) = fread(fid,1,'double');
    %     d=fread(fid,imSize,'ubit8');
    %     dd=uint8(fread(fid,imSize,'ubit8'));
    %     dd=uint8(fread(fid,imSize,'ubit8'));
    %     data=dd;
    %data = uint8(fread(fid,imSize,'ubit8'));
    frame = uint8(fread(fid,imSize,'uint8'));
    Data = reshape(frame, [h,w]);
    Data = flip(transpose(Data));
%     Data = flip(transpose(Data));
%     Data = flip(transpose(Data));
    Data=Data(floor(min(yi2)):ceil(max(yi2)), floor(min(xi2)):ceil(max(xi2)));
%     data=fread(fid,imSize,'uint8');
%     data=fread(fid,imSize,'uint8');
%     data=double(fread(fid,imSize,'ubit8'));
    im = Data;
    %im=reshape(Data,[w,h])';
    
    im2 = histeq(im); % increase contrast between pupil and iris
    im2=imadjust(im);
    %im2 = im;
    bwim=im2bw(im2,conthres); % convert to binary black/white using middle gray as the threshold
%         imshow(bwim)
    % CLIP
%     bwim=bwim(1:75,:);
    
    % 2nd variable here determines the range of radii to look for
    [center,radius] = imfindcircles(bwim,[10 20],'ObjectPolarity','Dark'); % change object polarity to bright/dark
    
%       imshow(im2)
%       viscircles(center, radius,'EdgeColor','b')
    
    %was using [.1*w .45*w] 
    if numel(radius)>1
        %hello
        % if multiple circles detected, use the one closest to the center
        % of the image
        dist_from_center = sqrt( (center(:,1) - w/2).^2 + (center(:,2) - h/2).^2);
        [~,min_idx] = min(dist_from_center);
        
        center = center(min_idx,:);
        radius = radius(min_idx,:);
        
%            hello
    end
    
    
    if(~isempty(radius))
        
        pupil_diam(iF) = 2*radius; %get diameter
        
        
        sizeChange = abs(prevDiam-pupil_diam(iF));
%         centerChange=abs(center(1)-prevCenter);
        
        
        if ctr>0 && sizeChange>diamCheck %if there is a sudden huge change in diameter, rerun imfindcircles with flipped polarity
            
            
%             hello
            bwim=im2bw(im2,conthres);
            [center,radius] = imfindcircles(bwim,[5 20],'ObjectPolarity','Dark'); % change object polarity to bright/dark
%            prevCenter=center(1);
            if numel(radius)>1
                
                % if multiple circles detected, use the one closest to the center
                % of the image
        dist_from_center = sqrt( (center(:,1) - w/2).^2 + (center(:,2) - h/2).^2);
        [~,min_idx] = min(dist_from_center);
        
      %  [~,min_idx]=min(radius);
        center = center(min_idx,:);
        radius = radius(min_idx,:);
            end
            newSizeChange=abs(prevDiam-2*radius);
            
            if(~isempty(radius))% && newSizeChange<diamCheck %if there is a circle and it's size is not 10 pixels off from previous
                pupil_diam(iF) = 2*radius;
            else
                pupil_diam(iF) = nan;
            end
            
        end
        
%         if isempty(radius)
%         hello
%         end
        
        if(showIm)
            figure(2); clf;
            imagesc(bwim) %was im
            colormap gray
            viscircles(center,radius,'color','g');
            text(0.75*h,0.2*w,sprintf('D = %.1f',2*radius),'FontSize',20,'Color','g');
            text(0.25*h,0.2*w,sprintf('t = %.1f',pupil_tstamp(iF)),'Fontsize',20,'Color','g');
            drawnow;
         
%             fim=figure(3);
%             subplot(131); imshow(im);
%             subplot(132); imshow(im2);
%             subplot(133); imshow(bwim);
%             pause
%             clear fim
            
            if(saveVideo)
                F(end+1) = getframe(gca);
            end
            
        end

    else
        pupil_diam(iF) = nan;
    end

    
    %use last measured pupil val to compare with
    if ~isnan(pupil_diam(iF))
        prevDiam=pupil_diam(iF);
    elseif isnan(pupil_diam(iF))
        prevDiam=prevDiam;
    end
    
    if ~isempty(radius)==1
    radiusall(iF)=radius;
    else
    radiusall(iF)=NaN;
    end
    
    ctr=ctr+1;
    iF/endloopwith
end
fclose(fid);

figure;
plot(pupil_tstamp,pupil_diam/2,'k-');
xlabel('Frames');
ylabel('Pixels');
title('Pupil Diameter');

if (savePupil)
save([date,'_Pupil','.mat'],'pupil_diam')
end

%%%%%%
if(saveVideo)
    F2 = F(2:end);
    v = VideoWriter([pupilVidName '.avi'],'Motion JPEG AVI');
    v.FrameRate = 30;
    open(v);
    writeVideo(v,F2);
    close(v);
end


%% synchronize pupil data
%%% NEED TO UPDATE, load in behavior data, and sync
% data = exper.M;
% %find synch pulse and throw out everything that occurred before that. This
% %puts neural and behav data in same start time
% 
% synchro=find(data(:,4)==1,2);
% synchro_tstamp = data(synchro,1);
% 
% 
% [~,pupil_start] = min(abs(pupil_tstamp-synchro_tstamp));
% %
% pupil_t = pupil_tstamp(pupil_start:end);
% pupil_d = pupil_diam(pupil_start:end);
% 
% pupil_t = pupil_t-(pupil_t(1));
% 
% [~,pupil_stop] = min(abs(pupil_t-endTime));
% pupil_d = pupil_d(1:pupil_stop);
% pupil_t = pupil_t(1:pupil_stop);
% 
% %%%% Note in this case since we're doing only 1 plane, imaging rate is
% %%%% actually higher than the ppupil vid framerate. So need to upsample.
% %%%% Normally this would be a downsample step.
% fs_up_pupil = length(pupil_t)/size(df,2);
% xq = linspace(1,length(pupil_d),length(t));
% pupil_resamp = interp1(1:length(pupil_d),pupil_d,xq);
% 
% % normalize pupil diam
% pupil = ( (pupil_resamp - min(pupil_resamp))/max(pupil_resamp - min(pupil_resamp)) )';
% 
% 
% 
% save(fullfile(currentFolder,outDir,'pupildiam.mat'),'pupil');
