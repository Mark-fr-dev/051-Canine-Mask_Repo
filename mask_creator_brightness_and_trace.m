clear

load d:\canine_study\051_Masks\bright_mask_BL_ele23_th30

mask_bright =mask(1:end-100,:);
load D:\Canine_study\051_Masks\00_BL0_Masks\BLmaskSAX_ele23_1
mask_hand =mask;

figure(21)
subplot(1,3,1)
imagesc(mask_hand)
title('hand')

subplot(1,3,2)
imagesc(mask_bright)
title('bright')

subplot(1,3,3)
imagesc(mask_bright.*mask_hand)
title('combine')
