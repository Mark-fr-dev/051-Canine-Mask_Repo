% Code: B-Mode brightness -> mask
% Input: B-mode axial x lateral matrix 
% Output: ROI mask based on B-Mode brightness
% Purpose: To produce a smooth 3D mask for K9 data display

clear


%% File Handling
run D:\Canine_study\Canine_study_file_names.m

%% Size of data (Canine in this case)
run D:\Canine_study\Canine_study_data_sizes.m



vid_file = sprintf('B_mode_brightness_masks_%s.avi', datestr(now, 'yyyy-mm-dd_HHMMSS'));
video_flag=0;
vv = open_vid(video_flag,vid_file);

file_dir= 'D:\Canine_study\051_Masks_repo\Output_data\';

disease_name = "HO"
ele_num = 20;
for frame=ED_frame:26
    for slice = BL_ele_slice_st:BL_ele_slice_st

        file_name = [file_dir 'HO_bmode_frame' num2str(frame) '_ele' num2str(ele_num)];
        load(file_name)

        [Na, Nl]=size(gray_data);

        for TH=30%20:50

        %Display B-Mode image
        figure(1);
        tiledlayout(1,2)
        nexttile
        imagesc(gray_data)
        colormap('gray')
        colorbar
        clim([0 255])
        title(strcat(disease_name, 'Bmode fr', num2str(frame), ' ele', num2str(ele_num)))

        %Manually selected ROI
        strat_a=120; end_a = 1458;
        strat_l=3; end_l = 60;
        
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
        title(strcat('Mask ', disease_name,  ' Bmode fr', num2str(frame), ...
            ' ele', num2str(ele_num), ' TH', num2str(TH)))

        capture_frame(video_flag,vv)
        savename = strcat("Output_data\",disease_name, " Bright mask fr ", num2str(frame),...
            " ele ", num2str(ele_num))
        save(savename,'mask');

        end
    end
end

close_vid(video_flag,vv)



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

function capture_frame(flag,vid_file)
    if(flag ==1)
        frame = getframe(gcf);
        writeVideo(vid_file,frame);
    else
        return
    end
end

function close_vid(flag, vid_file)
    if(flag ==1)
        close(vid_file)
    else
        return
    end
end

function video = open_vid(flag,file_name)
    if(flag==1)
        video = VideoWriter(file_name);
        video.FrameRate=2;
        open(video)
    else
        video = 0;
        return
    end
end
