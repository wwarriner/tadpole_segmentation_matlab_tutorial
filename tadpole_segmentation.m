%% NOTES

% In MATLAB, right-click a function name -> help, to see more info on it.

% If your images are larger than 320x240 pixels, the numeric values in ALL_CAPS
% may need to increase proportionally. If you double the image size, double the
% radius values, quadruple the area values.

% Try running this script on some images you have hand-calculated to make sure
% it lives up to your expectations. Bear in mind the units of area are pixel^2.
% You'll need to convert to metric based on your scale bar.

% Going forward, I recommend making sure the background has a high contrast with
% the foreground (yellow on blue is good!). I also recommend making sure there
% is space between every tadpole, at least enough to see the blue background
% between them.

% Once you understand what's going on in the script, I have an idea for how to
% estimate length and width (I've done something similar before in my PhD work
% for metal castings). I'll leave it to you to implement it, but we can discuss
% if you have questions. Some ideas to get you started:
%  - https://www.mathworks.com/help/images/ref/bwdist.html
%  - https://www.mathworks.com/help/images/ref/bwskel.html

%% SETUP
clear; % deletes variables in workspace
imtool close all; % closes all imtool windows
close all; % closes all figure windows

SHOW_IMAGES = true; % change to true to visualize intermediate steps

%% PREPARE IMAGE
file = fullfile("images", "1.tiff");
im_orig = imread(file);

% convert from [0, 255] range of uint8 to [-inf, +inf] range of double-precision
% floating point, enables some math later
im = double(im_orig) ./ 255;
if SHOW_IMAGES; imtool(im); end

%% TRANSFORM TO GRAYSCALE
red = rescale(im(:, :, 1));
if SHOW_IMAGES; imtool(red); end
green = rescale(im(:, :, 2));
if SHOW_IMAGES; imtool(green); end
blue = rescale(im(:, :, 3));
if SHOW_IMAGES; imtool(blue); end

% I happened to notice the background is blueish, the tadpoles are yellowish
% If we take their difference we should get high contrast between tadpoles
yellow = (red + green) / 2;
intensity = yellow - blue;
gray = rescale(intensity);
if SHOW_IMAGES; imtool(gray); end

%% ADJUST IMAGE
% top-hat transform tends to drop-out large blobs of high intensity, and keep
% small ones, this basically removes uneven background lighting if done right.
% Large means larger than the radius below, so we want something close to the
% half the width of the tadpoles.
RADIUS = 12;
gray = imtophat(gray, strel('disk', RADIUS));
gray = rescale(gray);
if SHOW_IMAGES; imtool(gray); end % compare with the original image!

% the histogram produced by this code shows log(# pixels) at each intensity
% value from 0 to 1. Most of the image is dark (the background, left side of
% histogram), but there is a hill of values in the brighter region to the right:
% those are the tadpoles.
if SHOW_IMAGES
    figure_handle = figure();
    axes_handle = axes(figure_handle);
    BIN_COUNT = 256; % number of gray levels in 8-bit image
    histogram(axes_handle, gray, BIN_COUNT);
    axes_handle.YScale = "log";
end

%% BINARIZE IMAGE
% graythresh() automatically picks a value for the threshold using image
% statistics (Otsu's method, specifically).
threshold_value = graythresh(gray);
binary_image = imbinarize(gray, threshold_value);
GREEN = [0 1 0];
if SHOW_IMAGES; imtool(labeloverlay(gray, binary_image, 'colormap', GREEN)); end

%% IMAGE MORPHOLOGY
% grow the tadpoles a bit to close up gaps
RADIUS = 1;
binary_image = imdilate(binary_image, strel('disk', RADIUS));
if SHOW_IMAGES; imtool(labeloverlay(gray, binary_image, 'colormap', GREEN)); end

% get rid of "small" regions that definitely aren't tadpoles
PIXEL_AREA = 25;
binary_image = bwareaopen(binary_image, PIXEL_AREA);
if SHOW_IMAGES; imtool(labeloverlay(gray, binary_image, 'colormap', GREEN)); end

% fill in holes within the tadpoles
binary_image = imfill(binary_image, 'holes');
if SHOW_IMAGES; imtool(labeloverlay(gray, binary_image, 'colormap', GREEN)); end

% shrink the tadpoles and remove "small" regions again
RADIUS = 5;
foreground_markers = imerode(binary_image, strel('disk', RADIUS));
PIXEL_AREA = 25;
foreground_markers = bwareaopen(foreground_markers, PIXEL_AREA);
if SHOW_IMAGES; imtool(labeloverlay(gray, foreground_markers, 'colormap', GREEN)); end

% skeletonize the background (hence the ~, it turns 1 to 0 and 0 to 1)
% also remove small pieces of the skeleton, to avoid segmenting curled tadpoles
% into multiple pieces
RADIUS = 3;
background_markers = imdilate(binary_image, strel('disk', RADIUS));
background_markers = bwskel(~background_markers);
background_markers = bwareaopen(background_markers, 25);
if SHOW_IMAGES; imtool(labeloverlay(gray, background_markers, 'colormap', GREEN)); end

% compute the image gradient, which is basically pixels where intensity changes
% suddenly: edges
gradient_magnitude = rescale(imgradient(gray));
if SHOW_IMAGES; imtool(imfuse(gradient_magnitude, gray)); end

% impose the markers onto the image to prepare for watershed
gradient_magnitude = imimposemin(gradient_magnitude, foreground_markers | background_markers);
gradient_magnitude(gradient_magnitude == -inf) = -1; % this is only needed to show the change visually
if SHOW_IMAGES; imtool(imfuse(gradient_magnitude, gray)); end

% pad the array with a gutter to avoid removing tadpoles with imclearborder()
% later on
PAD_AMT = 1; % no need to change me!
gradient_magnitude = padarray(gradient_magnitude, [PAD_AMT PAD_AMT], "replicate", "both");

% perform the watershed to recover segments
% this is challenging to explain in comments
% see: https://imagej.net/Classic_Watershed
ws = watershed(gradient_magnitude);

% clear the border to remove the background segment
ws = imclearborder(ws);

% remove the padding from earlier
ws = ws(PAD_AMT + 1 : end - PAD_AMT, PAD_AMT + 1 : end - PAD_AMT);

% relabeling to start at 1 and keep labels contiguous
% we don't want [2, 3, 5, 10] for example, should be [1, 2, 3, 4]
labels = unique(ws(ws > 0));
new_label = 1;
for label = labels.'
    ws(ws == label) = new_label;
    new_label = new_label + 1;
end

% dilate each segment to include the edges of the tadpoles
for label = labels.'
    expanded = ws == label;
    expanded = imdilate(expanded, conndef(2, 'minimal'));
    ws(expanded) = label;
end

% clear small segments that aren't part of tadpoles
PIXEL_AREA = 50;
retain = bwareaopen(ws, PIXEL_AREA);
ws(~retain) = 0;

%% OUTPUT
imtool(labeloverlay(im_orig, ws));

stats = regionprops(ws, "area");
stats = struct2table(stats);
summary_stats = grpstats(stats, "", ["min" "mean" "max"]);

fprintf("VALUES IN PIXEL UNITS:" + newline);
disp(summary_stats);
