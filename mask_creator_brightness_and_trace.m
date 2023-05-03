% Code: Brightness mask + hand mask -> composite mask
% Input: Brightness mask from K9_brightness_mask_for_3D and a handcut mask
% Output: Composite mask
% Purpose: To produce a smooth 3D mask for K9 data display
% Authors: Poshun and Mark de Villiers
% date: 20203/4/28

clear

%% File Handling
run D:\Canine_study\Canine_study_file_names.m
%% Size of data (Canine in this case)
run D:\Canine_study\Canine_study_data_sizes.m
%% Reusable functions
addpath("D:\Canine_study")


vid_file = sprintf('Composite_masks_%s.avi', datestr(now, 'yyyy-mm-dd_HHMMSS'));
video_flag=1;
vv = open_vid(video_flag,vid_file);

brt_file_dir= 'D:\Canine_study\051_Masks_repo\Output_data\';

for dis = 2%BL_disease:HO_disease
    if dis == BL_disease
        ele_st = BL_ele_slice_st;
        ele_end = BL_ele_slice_end;
        disease_name = 'BL';
        out_file_stub = strcat("Output_data\",disease_name, " Bright mask fr ");
        brt_msk_stub = strcat(brt_file_dir, disease_name ,'_Bright_mask_fr_');
        hand_msk_stub = strcat('D:\Canine_study\051_Masks\00_BL0_Masks\BLmaskSAX_ele');
    elseif dis ==  HO_disease
        ele_st = HO_ele_slice_st;
        ele_end = 23;
        disease_name = 'HO';
        out_file_stub = strcat("Output_data\",disease_name, " Bright mask fr ");
        brt_msk_stub = strcat(brt_file_dir, disease_name ,'_Bright_mask_fr_');
        hand_msk_stub = strcat('D:\Canine_study\051_Masks\01_HL0_Masks\HOmaskSAX_ele');
    end

    for slice=ele_st:ele_end
        for frame=ED_frame

            %load d:\canine_study\051_Masks\bright_mask_BL_ele23_th30
            file_name = [brt_msk_stub num2str(frame) '_ele_' num2str(slice)];
            load (file_name);

            if dis == 2
                mask_bright =mask(1:end-100,:);
            else
                mask_bright = mask;
            end

            file_name = [hand_msk_stub num2str(slice) '_1'];
            load(file_name);
            %load D:\Canine_study\051_Masks\00_BL0_Masks\BLmaskSAX_ele23_1
            mask_hand =mask;


            figure(21)
            subplot(1,3,1)
            imagesc(mask_hand)
            title('hand')

            subplot(1,3,2)
            imagesc(mask_bright)
            title('bright')


            mask_comb = mask_bright.*mask_hand;
            subplot(1,3,3)
            imagesc(mask_comb)
            title('combined')

            capture_frame(video_flag,vv)

            savename = strcat("Output_data\",disease_name, "_combined_mask_fr_", num2str(frame),...
                "_ele_", num2str(slice));
            save(savename,'mask_comb');
        end
    end
end
close_vid(video_flag,vv)



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