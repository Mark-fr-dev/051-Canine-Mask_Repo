%%%%% Data read for raw data
%%%%% and compare to the measurement result with the short axis view (at a
%%%%% certain axial depth).
%%%%% Created by Geng-Shi Jeng 7/2/2017
% 0728 provide to handwriteROI
%line 27 32 68 69 168need to fix
%% Modified to run on Mark de Villiers PC 3/8/2023


clear all
addpath D:\Canine_Data\code_STMat
%% File Handling
run D:\Canine_study\Canine_study_file_names.m

%% Size of data (Canine in this case)
run D:\Canine_study\Canine_study_data_sizes.m

intp_factor = 1;  % interpolation factor in lateral and elevational directions

% rfDSEA16HO0
%seq_folder= 'Canine';
%seq_folder_sub = 'originaldata';
disease_name = BL_disease_name;
%disease_name = HO_disease_name;
base_name ='rfDSEA16HO0'
seq_id = 'rfDSEA16HO0';   % The folder name and the filename
N_file_index = 31;         % Number of files in the fholder
%rf_folder = ['D:\Canine_Data\canine_data\' seq_id   '\'  seq_folder_sub];

%%%% Extract one of the raw data to know the imaging parameters
% seq_folder_sub = 'originaldata';
% rf_folder = ['F:\Research\SpeckleTracking\Canine\' seq_folder '\' base_name '\' seq_folder_sub];
rf_fname = ['D:\Canine_Data\canine_data\rfDSEA16HO0\originaldata\' base_name num2str(1) '.mat'];
load(rf_fname)

c           = 1.54;     % [mm/us]
r           = s.Time_usec  * c / 2;
dr          = r(2) - r(1);  % mm
ThetaDegs   = fliplr(s.LatAngle_degs);  % lateral spanning degrees in polar coordinate
ThetaRads	= ThetaDegs * pi/180;
ThetaRads   = interp(ThetaRads,intp_factor); 
dthd        = min(diff(ThetaDegs));

PhiDegs     = fliplr(s.EleAngle_degs);  % elevational spanning degrees in polar coordinate
PhiRads     = PhiDegs * pi/180;
PhiRads     = interp(PhiRads,intp_factor);
dphid       = min(diff(PhiDegs));



%%% Segmentation and compare to ground truth data

   filt_r=5;%3				;%Filter limit (rows)
   filt_c=12			    ;%Filter limit (cols)
   filt_e=5;%3				;%Filter limit (elevation)
   

% load mask_f1
% 
% mask_tmp = mask3(filt_c+1:end-filt_c,filt_r+1:end-filt_r);
% 
% mask3 = mask_tmp;
% 
% clear mask_tmp
% 

% ele_slic is ROI ele slice 
%ele_slice = 1;
if disease_name == HO_disease_name
    axi_s = 500; axi_end = 2360;
     Msk_dir = Msk_data_dir_HO;


elseif disease_name == BL_disease_name
    axi_s = 600; axi_end = 2360;
     Msk_dir = Msk_data_dir_BL;
end



for ele_slice = 8:48
    for ff = 1
    
  rfdata_dir = ['D:\Canine_Data\canine_data\' seq_id  '_sigbb'];
  %%% load pre-deformation data 
  fnamesave = sprintf('%s\\%s_truc_nointp__%d',rfdata_dir,base_name,ff);
  load(fnamesave);
  PreImag = rfdata;

  %%% load post-deformation data
  fnamesave = sprintf('%s\\%s_truc_nointp__%d',rfdata_dir,base_name,ff+1);
  load(fnamesave);
  DefImag = rfdata;

  [axi_N,lat_N,ele_N] = size(DefImag);

  clear rfdata
  
  if ff == ff
     MAX_RF =  max(max(max(abs(PreImag))));
  end
         
       
           azi_factor = intp_factor;
           data_temp2 = zeros(axi_N,lat_N*azi_factor,ele_N);
        
           for gg=1:ele_N
             %for ll=1:axi_N
               %data_temp2(ll,:,gg) = double(interp(squeeze(sig_bb(ll,:,gg)),azi_factor));
             data_temp2(:,:,gg)=double(interpft(squeeze(PreImag(:,:,gg)),lat_N*azi_factor,2));
             %end
           end

           zen_factor = 1;
           [axi_N,lat_N,ele_N] = size(data_temp2);
           data = zeros(axi_N,lat_N,ele_N*zen_factor);
           
           for gg=1:lat_N
             %for ll=1:axi_N
              %data(ll,gg,:)=double(interp(squeeze(data_temp2(ll,gg,:)),zen_factor));
              data(:,gg,:)=double(interpft(squeeze(data_temp2(:,gg,:)),ele_N*zen_factor,2));
             %end
           end
           
           PreImag = data;
           
           clear data_temp2 data 
  
   [axi_N, lat_N, ele_N] = size(PreImag);
  %figure
  DR = 50;
  
  lpf = [1 3 3 1];
  lpf = lpf/sum(lpf);
  
  logdata = 20*log10(abs(PreImag(:,:,ele_slice))/MAX_RF);
  gray_data = (logdata+DR)/DR*255;
  
  
%  file_dir = ['E:\' seq_id '_0628_2ptf\']
% file_name= [file_dir   seq_id  '_0628_2ptf_result_disp_' num2str(ff) num2str(ff+1)  ];
% load(file_name)
% axial_disp = squeeze(axial_disp(:,:,ele_slice));
% lat_disp = squeeze(lat_disp(:,:,ele_slice));
% ele_disp = squeeze(ele_disp(:,:,ele_slice));
% NCC = squeeze(NCC(:,:,ele_slice));
% 
% shift_ry_org_seg_med = medfilt2(lat_disp,[20 5]);
% shift_ez_org_seg_med = medfilt2(ele_disp,[20 5]);
% shift_cx_org_seg_med = medfilt2(axial_disp,[20 5]);
% max_v_intp_org = medfilt2(NCC,[20 5]);

  
         %%%log compression
        figure(9)
        
%         imagesc(NCC)  
%         xlabel('Lateral (Degrees)')
%         ylabel('Axial (mm)')
%         caxis([0 1])
%         colorbar
%         title(['Frame #' num2str(ff)])

        image(ThetaDegs,r(axi_s:end),gray_data)  
        xlabel('Lateral (Degrees)')
        ylabel('Axial (mm)')
        colormap(gray(256))
        title(['Frame #' num2str(ff) 'Slice #' num2str(ele_slice)])

  
  gray_data = (gray_data > 0).*gray_data;
  gray_data_ds = gray_data ;
  imagename = ['Output_data\' seq_id '_frame' num2str(ff) '_ele' num2str(ele_slice) '.png'];
  imwrite(gray_data_ds , gray(256),imagename)  
  
  save_name = ['Output_data\' disease_name '_bmode' '_frame' num2str(ff) '_ele' num2str(ele_slice)];
  save(save_name,'gray_data')

  
  %% Segmentation 
%   mask1 = roipoly;
%   mask2 = roipoly;
%   mask3 = abs(mask1 - mask2);
   % mask3 = roipoly;
  
  
%   mask = mask3(1:axial_range,1:20:end);
%   centroid = regionprops(true(size(mask)), mask,  'WeightedCentroid');
%   
%   %save file 
%   savefilename = [ Msk_dir disease_name 'maskSAX_ele' num2str(ele_slice) '_' num2str(ff) '.mat']
%   save(savefilename,'mask', 'centroid') 
  
%   mask_tmp = mask3(filt_c+1:end-filt_c,filt_r+1:end-filt_r);
%   mask3 = mask_tmp;
%   clear mask_tmp
  
  
%   mask = zeros(size(gray_data));
%   mask(1:end-1,1:end-1) = 1;
%   seg_data = activecontour(gray_data,mask,300);
%   
%   figure
%   image(seg_data)
%   colormap(gray(2))
  
  
  
  %%% Load displacement estimation result 

%   fnamesave2 = [base_name '_PMSTTF_d4_181_535_1_' num2str(ff)]; 
%    
%   load(fnamesave2);
% 
%   %mask3 = ones(size(shift_ry));
%   shift_ry_seg_med = medfilt2(shift_ry.*mask3,[20 5])*dthd;
%   shift_ez_seg_med = medfilt2(shift_ez.*mask3,[20 5])*dphid;
%   shift_cx_seg_med = medfilt2(shift_cx.*mask3,[20 5])*dr;
% 
% 
%   shift_ry_abs_seg_med = medfilt2(shift_ry_abs.*mask3,[20 5])*dthd;
%   shift_ez_abs_seg_med = medfilt2(shift_ez_abs.*mask3,[20 5])*dphid;
%   shift_cx_abs_seg_med = medfilt2(shift_cx_abs.*mask3,[20 5])*dr;
% 
% 
%   shift_ry_ntf_seg_med = medfilt2(shift_ry_ntf.*mask3,[20 5])*dthd;
%   shift_ez_ntf_seg_med = medfilt2(shift_ez_ntf.*mask3,[20 5])*dphid;
%   shift_cx_ntf_seg_med = medfilt2(shift_cx_ntf.*mask3,[20 5])*dr;
%   
%   
%   %%% Load displacement with original Linux codes (with correlation filter
%   %%% only)
%   
%   
%   fname = ['F:\Research\SpeckleTracking\Canine\DSEA16\' base_name '\output\' base_name sprintf('_f%03d-%03d_s%03d.out',ff,ff+1,(56/2-1))]; 
%   
%   if (exist('fname') ~=1),
%   fname = input('Filename : ','s');	% Input filename
% end
% 
% if (exist(fname) ~=2),
%   disp(['File ' fname ' does not exist'])
%   %clear fname
%   return
% end
% 
% disp(['Reading file ' fname ' ...']);
% 
% 
% % Check if valid Binary file
% %fid = fopen(fname,'r');  % accumulation
% endian = 'b';
% fid = fopen(fname,'r',endian);
% 
% head=fread(fid,8,'char');
% head = setstr(head');
% valid_bin_old = (strcmp (head,'#CORR v1'));
% valid_bin = (strcmp (head,'DISP01.0'));
% valid_mat = 0;
% if ((~valid_bin_old)&(~valid_bin_old)&(length(fname)>4))
%   valid_mat=(strcmp (fname((length(fname)-3):length(fname)),'.mat'));
% end
% 
% 
% if (valid_bin)
%   win = fscanf(fid,'%d',6);
%   win_x = win(1):win(2):win(3);
%   win_y = win(4):win(5):win(6);
%   Nx = length(win_x);
%   Ny = length(win_y);
%   
%   % find beginning of data
%   data=fgets(fid);
%   currpos=ftell(fid);      %add
%   while (data~=(-1) & ~(any(data==12))),
%     oldpos = currpos;      %add
%     data=fgets(fid)
%     currpos = ftell(fid)  %add
%   end
% 
%   % eyw -- occasionally, file position indicator moves an extra step
%   % (3 steps instead of 2) on the "\f\n" delimiting line, so move it back
%   step=currpos-oldpos;
%   ftell(fid);
%   if(step~=2)        % size(data,2) is sometimes 2, sometimes 3
%     stepdiff=2-step;
%     fseek(fid,stepdiff,'cof')   % move current position
%   end
%   ftell(fid);
%   %--  
%   
%   % Read data
%   [data,count]=fread(fid,[4*Nx,Ny],'float');
%   fclose(fid);
% 
%   
%   if (count ~= (Nx*Ny*4)) 
%     disp(sprintf('Warning 1: Only %d out of %d values read',count/4,Nx*Ny))
%   end
%   
%   clear fname;
%   
%   % Reshape vectors into matricies
%   disp('Reshaping matricies...');
%   R = data(1:4:(4*Nx),:).';
%   x_shift = data(2:4:(4*Nx),:).';
%   y_shift = data(3:4:(4*Nx),:).';
%   z_shift = data(4:4:(4*Nx),:).';
% 
%   %clear data count win fid head valid_bin_old valid_bin valid_mat Nx Ny
% 
%   
% elseif (valid_bin_old)
%   win = fread(fid,6,'short');
%   win_x = win(1):win(2):win(3);
%   win_y = win(4):win(5):win(6);
%   Nx = length(win_x);
%   Ny = length(win_y);
%   [data,count]=fread(fid,[4*Nx,Ny],'float');
%   fclose(fid);
% 
%   if (count ~= (Nx*Ny*4)) 
%     disp(sprintf('Warning 2: Only %d out of %d values read',count/4,Nx*Ny))
%   end
%   
%   clear fname;
% 
%  
%   
%   % Reshape vectors into matricies
%   disp('Reshaping matricies...');
%   R = data(1:4:(4*Nx),:)';
%   x_shift = data(2:4:(4*Nx),:)';
%   y_shift = data(3:4:(4*Nx),:)';
%   z_shift = data(4:4:(4*Nx),:)';
% 
%   %clear data count win fid head valid_bin_old valid_bin valid_mat Nx Ny
%   
% elseif (valid_mat)
%   fclose(fid);
%   clear R
%   eval(['load ' fname]);
% 
%   clear fname;
%   Nx = length(win_x);
%   Ny = length(win_y);
% 
%   % Reshape vectors into matricies
%   disp('Reshaping matricies...');
%   if (exist('R') == 1)
%     [tmp_r,tmp_c] = size(R);
%     if ((tmp_r~=Ny) & (tmp_c~=Nx))
%       R = reshape(R,Nx,Ny)';
%     end
%   end
%   [tmp_r,tmp_c] = size(x_shift);
%   if ((tmp_r~=Ny) & (tmp_c~=Nx))
%     x_shift = reshape(x_shift,Nx,Ny)';
%   end
%   [tmp_r,tmp_c] = size(y_shift);
%   if ((tmp_r~=Ny) & (tmp_c~=Nx))
%     y_shift = reshape(y_shift,Nx,Ny)';
%   end
%  
%   %clear tmp_r tmp_c fid head valid_bin_old valid_bin valid_mat Nx Ny
%   
% else
%   keyboard
%   fclose(fid);
%   disp('Not a valid data file');
%   %clear fname fid head valid_bin_old valid_bin valid_mat
%   return
% end  
% 
% 
% y_shift = y_shift.';
% z_shift = z_shift.';
% x_shift = x_shift.';
% R = R.';
% 
% shift_ry_org = y_shift(700-301+1+filt_c:2000-301+1-filt_c, 9+filt_r:55-filt_r)*dthd;
% shift_ez_org = z_shift(700-301+1+filt_c:2000-301+1-filt_c, 9+filt_r:55-filt_r)*dphid;
% shift_cx_org = x_shift(700-301+1+filt_c:2000-301+1-filt_c, 9+filt_r:55-filt_r)*dr;
% max_v_intp_org = R(700-301+1+filt_c:2000-301+1-filt_c, 9+filt_r:55-filt_r);
% 
% shift_ry_org_seg_med =shift_ry_org_seg_med.*mask3;
% shift_ez_org_seg_med = shift_ez_org_seg_med.*mask3;
% shift_cx_org_seg_med = shift_cx_org_seg_med.*mask3;
% max_v_intp_org =max_v_intp_org.*mask3;
% 
% 
% %%% Regular speckle tracking with correlation filter only 
% figure
% subplot(2,2,1)
% imagesc(ThetaDegs,r(600:end) , shift_ry_org_seg_med)
% %colormap(jet)
% colorbar
% title('Lateral')
% xlabel('Lateral(degrees)')
% ylabel('Axial(mm)')
% % caxis([-1*dthd 1*dthd])
% caxis([-2 2])
% 
% subplot(2,2,2)
% imagesc(ThetaDegs,r(600:end) ,shift_ez_org_seg_med)
% colorbar
% title('Elevational')
% xlabel('Lateral(degrees)')
% ylabel('Axial(mm)')
% % caxis([-1*dphid 1*dphid])
% caxis([-2 2])
% 
% subplot(2,2,3)
% imagesc(ThetaDegs,r(600:end) ,shift_cx_org_seg_med)
% colorbar
% title('Axial')
% xlabel('Lateral(degrees)')
% ylabel('Axial(mm)')
% % caxis([-5*dr 5*dr])
% caxis([-10 10])
% 
%  subplot(2,2,4)
% imagesc(ThetaDegs,r(600:end) ,max_v_intp_org)
% colorbar
% caxis([0 1])
% xlabel('Lateral(degrees)')
% ylabel('Axial(mm)')
% title('NCC')
% 
% % 
% % %%% PM and with correlation filter only 
% % figure
% % subplot(2,2,1)
% % imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_ry_ntf_seg_med)
% % colorbar
% % title('Lateral')
% % xlabel( 'Lateral(degrees)')
% % ylabel('Axial(mm)')
% % caxis([-1*dthd 1*dthd])
% % 
% % subplot(2,2,2)
% % imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_ez_ntf_seg_med)
% % colorbar
% % title('Elevational')
% % xlabel( 'Lateral(degrees)')
% % ylabel('Axial(mm)')
% % caxis([-1*dphid 1*dphid])
% % 
% % 
% % subplot(2,2,3)
% % imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_cx_ntf_seg_med)
% colorbar
% title('Axial')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-5*dr 5*dr])
% 
% ax4 = subplot(2,2,4)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,max_v_intp_ntf)
% colormap(ax4,hot)
% colorbar
% caxis([0 1])
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% title('NCC')
% 
% 
% %%% PM with abs(NCC)
% figure
% subplot(2,2,1)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_ry_abs_seg_med)
% %colormap(jet)
% colorbar
% title('Lateral')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-1*dthd 1*dthd])
% 
% 
% subplot(2,2,2)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_ez_abs_seg_med)
% colorbar
% title('Elevational')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-1*dphid 1*dphid])
% 
% subplot(2,2,3)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_cx_abs_seg_med)
% colorbar
% title('Axial')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-5*dr 5*dr])
% 
% ax5 = subplot(2,2,4)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,max_v_intp_abs)
% colormap(ax5,hot)
% colorbar
% caxis([0 1])
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% title('NCC')
% 
% 
% %%% PM with tilt filtering
% figure
% subplot(2,2,1)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_ry_seg_med)
% colorbar
% title('Lateral')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-1*dthd 1*dthd])
% 
% 
% subplot(2,2,2)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_ez_seg_med)
% colorbar
% title('Elevational')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-1*dphid 1*dphid])
% 
% subplot(2,2,3)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,shift_cx_seg_med)
% colorbar
% title('Axial')
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% caxis([-5*dr 5*dr])
% 
% ax6 = subplot(2,2,4)
% imagesc(ThetaDegs(9+filt_r:55-filt_r),r(700+filt_c:2000-filt_c) ,max_v_intp)
% colormap(ax6,hot)
% colorbar
% caxis([0 1])
% xlabel( 'Lateral(degrees)')
% ylabel('Axial(mm)')
% title('NCC')
% 



end
end









           
           
           





