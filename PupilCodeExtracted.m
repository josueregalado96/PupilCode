%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%set up camera! for Initiationalization %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if vr.trackpupilday==1
% % %PUPIL CODE setup start:
NET.addAssembly('C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
vr.pupildirectory='C:\Users\RajasethupathyLab\Desktop\Pupil Data\';
vr.pupildirectory= [vr.pupildirectory datestr(now,'mmddyyyy') '\'];
vr.pupildirectoryimages=[vr.pupildirectory, datestr(now,'HHMMSS') '\'];

mkdir(vr.pupildirectory)
% mkdir(vr.pupildirectoryimages)
vr.pupilfilename = [vr.pupildirectory, strcat(vr.name, '_Pupil_',vr.filename), '.bin']

% Camera object handle
vr.cam = uc480.Camera;
% Open camera
vr.cam.Init(0);
vr.cam.Display.Mode.Set(uc480.Defines.DisplayMode.DiB);
vr.cam.Trigger.Set(uc480.Defines.TriggerMode.Software);
vr.cam.PixelFormat.Set(uc480.Defines.ColorMode.Mono8);

% %Set PixelClock to Max
% [~,vr.Min, vr.Max, ~]=vr.cam.Timing.PixelClock.GetRange()
% vr.cam.Timing.PixelClock.Set(vr.Min);
% %Set FrameRate to Min
% [~,vr.Min, vr.Max, ~]=vr.cam.Timing.Framerate.GetFrameRateRange()
% vr.cam.Timing.Framerate.Set(vr.Min);
% % idk
[~,vr.MemId] = vr.cam.Memory.Allocate(true);
[~,vr.Width,vr.Height,vr.Bits,~] = vr.cam.Memory.Inquire(vr.MemId);
% Set Exposure to Max
[~,vr.Min, vr.Max, ~] =vr.cam.Timing.Exposure.GetRange();
vr.cam.Timing.Exposure.Set(vr.Max)

vr.cam.Acquisition.Freeze(uc480.Defines.DeviceParameter.Wait)
[~,vr.tmp] = vr.cam.Memory.CopyToArray(vr.MemId);
vr.Data = reshape(uint8(vr.tmp),[vr.Width,vr.Height]);
vr.Data=imrotate(vr.Data,-90);
vr.Data=flipdim(vr.Data,2);
vr.Data = squeeze(vr.Data);
imwrite(vr.Data,'firstimg.tiff');
vr.roipic=imread('firstimg.tiff');
[vr.ROI,vr.x,vr.y]=roipoly(vr.roipic);
vr.x1=floor(min(vr.x));
vr.x2=ceil(max(vr.x));
vr.y1=floor(min(vr.y));
vr.y2=ceil(max(vr.y));

%set up AOI
vr.x1=roundn(vr.x1,1)-10;
vr.x2=roundn(vr.x2,1)+10;
vr.y1=roundn(vr.y1,1)-10;
vr.y2=roundn(vr.y2,1)+10;
[vr.x1 vr.x2 vr.y1 vr.y2]
vr.width=vr.x2-vr.x1;
vr.height=vr.y2-vr.y1;
if rem(vr.width,4)==2
    vr.width=vr.width+10;
end

[a,b,c]=vr.cam.Size.AOI.GetSizeRange;
[a,lower_left_x,lower_left_y,w,h] = vr.cam.Size.AOI.Get;
% cam.Size.AOI.Set(x1,y1,1000,1000)
vr.cam.Size.AOI.Set(vr.x1,vr.y1,vr.width,vr.height);
[a,lower_left_x,lower_left_y,w,h] = vr.cam.Size.AOI.Get;

%Set PixelClock to Max
[a,MinPC, MaxPC, Inc]=vr.cam.Timing.PixelClock.GetRange();
% vr.cam.Timing.PixelClock.Set(5);
%Set FrameRate to Min
[a,MinFR, MaxFR, Inc]=vr.cam.Timing.Framerate.GetFrameRateRange();
vr.cam.Timing.Framerate.Set(23);
% Set Exposure to Max
[a,MinE, MaxE, Inc] =vr.cam.Timing.Exposure.GetRange();
vr.cam.Timing.Exposure.Set(44);

[a,vr.bFR]=vr.cam.Timing.Framerate.Get();
[a,vr.bPC]=vr.cam.Timing.PixelClock.Get();
[a,vr.bE]=vr.cam.Timing.Exposure.Get();

vr.fid_pupil= fopen(vr.pupilfilename,'w');
fwrite(vr.fid_pupil,vr.width,'double');
fwrite(vr.fid_pupil,vr.height,'double');
fwrite(vr.fid_pupil,vr.Bits,'double');

vr.cam.Acquisition.Capture(uc480.Defines.DeviceParameter.DontWait);
vr.countcapture=0;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%then for each interation of virmen we start capturing images on the
%%%%50th iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if vr.iterations==50 || vr.iterations==5000 || vr.iterations==10000
    vr.syncOut = 1;
    %I think elements in the second argument are values to be sent on each of
    %the lines defined above
    vr.M(vr.iterations,3)=1;
    %PULSE SHOULD HAVE A CORRESPONDING VALUE SAVED IN VIRMEN
    %MAYBE A VR.PULSETIME VARIABLE, OUTPUT ON EVERY TELEPORT?
    if vr.trackpupilday==1
    vr.startrecordingpupil=1;
    end
end

if vr.startrecordingpupil==1
    % [~,vr.MemId] = vr.cam.Memory.Allocate(true);
    % [~,vr.Width,vr.Height,vr.Bits,~] = vr.cam.Memory.Inquire(vr.MemId);
    vr.countcapture = vr.countcapture+1;
    [~,vr.tmp] = vr.cam.Memory.CopyToArray(vr.MemId);
    % reformat image
    vr.Data = reshape(uint8(vr.tmp),[vr.Width,vr.Height]);
    vr.Data=imrotate(vr.Data,-90);
    vr.Data=flipdim(vr.Data,2);
    vr.Data = squeeze(vr.Data);
    vr.ROIy = 1:vr.height;
    vr.ROIx = 1:vr.width;
    vr.Data = vr.Data(vr.ROIy,vr.ROIx);
    fwrite(vr.fid_pupil,vr.countcapture,'double');
    fwrite(vr.fid_pupil,vr.Data(:));
    
%imwrite(vr.Data,[vr.pupildirectoryimages, sprintf('%d.png',vr.countcapture)])

disp(strcat('---------------------------------------------------------framerate===',num2str(vr.bFR)));
disp(strcat('---------------------------------------------------------expo===',num2str(vr.bE)));
disp(strcat('---------------------------------------------------------pixelcount===',num2str(vr.bPC)));
%vr.pupilarray(:,:,vr.countcapture)=vr.Data;
% display image
%  imshow(imadjust(vr.Data));
%  colormap gray
%  drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%then to finish the code we run this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
if vr.trackpupilday==1
fclose(vr.fid_pupil)
vr.cam.Acquisition.Stop(uc480.Defines.DeviceParameter.DontWait);
vr.cam.Exit;
end
