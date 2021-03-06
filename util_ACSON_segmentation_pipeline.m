function util_HM_SBEM_segmentation_pipeline(opt)


r_im = opt.read_im;
s_address = opt.save_address;
r_mask = opt.read_mask;
myelin_inx = opt.myelin_inx;
similarity_th = opt.similarity_th;
min_max_myelin_volume = opt.min_max_myelin_volume;
min_max_lbl_volume = opt.min_max_lbl_volume;
small_volume_th = opt.closing_small_volume_threshold;
bg = opt.background;
fg = opt.foreground;


%% Load

% Load stack of aligned images. Raw image should be between [0 255]
t = load(r_im); fields = fieldnames(t); raw_im = t.(fields{1}); raw_im = squeeze(raw_im);

% After alignment, if the output was not in the same size as the raw image stack, 
% a mask should define non-tissue voxels. 
% The mask is a binary 3D image with the same size as 
% the aligned image stack, equal to one for the non-tissue voxels and zero elsewhere. 
% After alignment, if the non-tissue voxels were cropped, mask would be empty. 


if isempty(r_mask)
    mask = [];
else
    t = load(r_mask);
    fields = fieldnames(t); mask = t.(fields{1});
end
clear t fields


%% BM4D filtering

if opt.filt_im == 0
    [filt_im,~] = bm4d(raw_im);
    save(strcat(s_address,'filt_im'),'filt_im','-v7.3')
    clear raw_im
else
    filt_im = raw_im;
    clear raw_im  
end


%% Canny edge detection
% STD Gaussian filter: sqrt(2).
% Weak and strong edges:0.25 and 0.6 times the maximum gradient magnitude.

filt_im = double(filt_im); filt_im(mask) = 255;
bound_im = false(size(filt_im));
for i = 1:size(filt_im,3)
    bound_im(:,:,i) = edge(filt_im(:,:,i),'canny');
end
bound_im(mask) = true;
bound_im = imdilate(bound_im,true(3));


%% Myelin detection

filt_im = padarray(filt_im/255, [1 1 1], 1);
bound_im = padarray(bound_im, [1 1 1], 1);
[r1,c1,h1] = size(filt_im);

myelin_inx = sub2ind([r1,c1,h1],myelin_inx(1),myelin_inx(2),myelin_inx(3));
myelin_rgn = false(r1,c1,h1);

if isinf(min_max_myelin_volume(2))
    min_max_myelin_volume(2) = numel(filt_im);
end

[myelin_rgn_inx,~] = util_regionGrowing(filt_im,myelin_inx,similarity_th,min_max_myelin_volume,bound_im,bg);
myelin_rgn(myelin_rgn_inx) = true;
            
myelin_rgn = imclose(myelin_rgn,true(3));
label = double(myelin_rgn);

bound_im(myelin_rgn) = true;

myelin_rgn = myelin_rgn(2:end-1,2:end-1,2:end-1);
save(strcat(s_address,'mat_myelin_rgn'),'myelin_rgn','-v7.3')

clear myelin_rgn


%% Closing small volumes

bw_close = bwareaopen(~bound_im,small_volume_th);
small_volumes = ~(bw_close|bound_im);
bound_im(small_volumes) = true;

small_volumes = small_volumes(2:end-1,2:end-1,2:end-1);
save(strcat(s_address,'mat_small_volumes'),'small_volumes','-v7.3')

clear small_volumes bw_close


%% Finding the location of seeds

rgnmaxima = zeros(r1,c1,h1);
for i = 2:h1-1
    D = bwdist(bound_im(:,:,i));
    rgnmaxima(:,:,i) = imregionalmax(D,4);
end
rps = find(rgnmaxima);
rps = rps'; 
clear rgnmaxima


%% regionGrowing

lbl = 1;
for seed = rps
    if bound_im(seed)==0
        [rgn_inx,flag] = util_regionGrowing(filt_im,seed,similarity_th,min_max_lbl_volume,bound_im,fg);
        if  flag==0
            label(rgn_inx) = lbl;
            bound_im(rgn_inx) = true;
            lbl = lbl+1;
        end
    end
end

label = label(2:end-1,2:end-1,2:end-1);
save(strcat(s_address,'label'),'label','-v7.3')


