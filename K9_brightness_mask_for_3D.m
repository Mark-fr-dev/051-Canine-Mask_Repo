% Code: B-Mode brightness -> mask
% Input: B-mode axial x lateral matrix 
% Output: ROI mask based on B-Mode brightness
% Purpose: To produce a smooth 3D mask for K9 data display

clear

frame_num = 1;
ele_num=23;

vid_file = sprintf('B_mode_brightness_masks_%s.avi', datestr(now, 'yyyy-mm-dd_HHMMSS'));
video_flag=1;
if(video_flag==1)
    vv= VideoWriter(vid_file);
    vv.FrameRate=2;
    open(vv)
end

file_dir= 'D:\Canine_study\051_Masks_repo\Output_data\';
file_name = [file_dir 'BL_bmode_frame' num2str(frame_num) '_ele' num2str(ele_num)];
load(file_name)

[Na, Nl]=size(gray_data);

for TH=20:50

%Display B-Mode image
figure(1);
tiledlayout(1,2)
nexttile
imagesc(gray_data)
colormap('gray')
colorbar
clim([0 255])
title(['BL Bmode fr' num2str(frame_num) ' ele' num2str(ele_num)])

%Manually selected ROI
strat_a=120; end_a = 1458;
strat_l=13; end_l = 51;

filtered_gray = medfilt2(gray_data,[34 7]);
mask =zeros(Na,Nl);
%TH = 35;  %B-Mode threshold brightness value 0-255
% TH_index = (gray_data>TH);
TH_index = (filtered_gray>TH);
mask(TH_index)=1;

%Ignore image outside this rectangular ROI
mask(1:strat_a,:)=0;mask(end_a:Na,:)=0;
mask(:,1:strat_l)=0;mask(:,end_l:Nl)=0;
mask = medfilt2(mask,[32 5]);


nexttile
imagesc(mask)
%colormap("parula")
colorbar
clim([0 1])
title(['Mask BL Bmode fr' num2str(frame_num) ' ele' num2str(ele_num) ' TH' ...
    num2str(TH)])


if(video_flag ==1)
        frame = getframe(gcf);
        writeVideo(vv,frame)
end

end


if(video_flag ==1)
    close(vv)
end

% Edge detection
% BW1 = edge(mask,'sobel');
% BW2 = edge(mask,'canny');
% 
% figure(2)
% tiledlayout(1,2)
% 
% nexttile
% imagesc(BW1)
% title('Sobel Filter')
% 
% nexttile
% imagesc(BW2)
% title('Canny Filter')


