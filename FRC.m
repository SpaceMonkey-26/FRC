clc; clear;
%% INPUTS
inputImageLocation1 = "Downloads\MB1.1_4a9p65z_x1N_OS-SIM_RL0_0o90o.tif";
inputImageLocation2 = "Downloads\MB1.1_4a9p65z_x1N_OS-SIM_RL0_45o135o.tif";
pixelSize = 0.065; % Pixel size: micrometers per pixel
binPercent = 0.02; % Percentage of total bins for smoothing: 0.02-0.05
window = 3; % Minimum window size for smoothing: 3, 5, or 7
%% IMPORT IMAGE DATA
imageData1 = imfinfo(inputImageLocation1);
imageData2 = imfinfo(inputImageLocation2);
zSlices = numel(imageData1);
height = imageData1(1).Height;
width  = imageData1(1).Width;
if (numel(imageData2) == zSlices) && (imageData2(1).Height == height) && (imageData2(1).Width == width)
    imageStack1 = zeros(height, width, zSlices);
    imageStack2 = zeros(height, width, zSlices);
    for z = 1:zSlices
        imageStack1(:,:,z) = double(imread(inputImageLocation1, z));
        imageStack2(:,:,z) = double(imread(inputImageLocation2, z));
    end
else
    error('Image #1 and Image #2 must have the same dimensions: x, y, & z');
end
%% COMPUTE FFT, PSD, & CC
imageFFT1 = zeros(height, width, zSlices);
imageFFT2 = zeros(height, width, zSlices);
imagePSD1 = zeros(height, width, zSlices);
imagePSD2 = zeros(height, width, zSlices);
imageCC = zeros(height, width, zSlices);
for z = 1:zSlices
    imageFFT1(:, :, z) = fftshift(fft2(imageStack1(:, :, z)));
    imageFFT2(:, :, z) = fftshift(fft2(imageStack2(:, :, z)));
    imagePSD1(:, :, z) = abs(imageFFT1(:, :, z)).^2;
    imagePSD2(:, :, z) = abs(imageFFT2(:, :, z)).^2;
    imageCC(:, :, z) = real(imageFFT1(:, :, z) .* conj(imageFFT2(:, :, z)));
end
%% COMPUTE CIRCULAR AVERAGES
maxRadius = (height + width) / 4;
binCount = round(3 * maxRadius / 4);
binCounter = zeros(binCount, 1, zSlices);
binValues = zeros(binCount, 3, zSlices);
yCenter = height / 2;
xCenter = width / 2;
for z = 1:zSlices
    for x = round(xCenter - maxRadius):round(xCenter + maxRadius)
        for y = round(yCenter - maxRadius):round(yCenter + maxRadius)
            if (x < 1) || (x > width) || (y < 1) || (y > height)
                continue
            end
            radius = sqrt((x - xCenter).^2 + (y - yCenter).^2);
            if radius > maxRadius
                continue
            end
            bin = min(floor((radius / maxRadius) * binCount) + 1, binCount);
            binValues(bin, 1, z) = binValues(bin, 1, z) + imagePSD1(y, x, z);
            binValues(bin, 2, z) = binValues(bin, 2, z) + imagePSD2(y, x, z);
            binValues(bin, 3, z) = binValues(bin, 3, z) + imageCC(y, x, z);
            binCounter(bin, 1, z) = binCounter(bin, 1, z) + 1;
        end
    end
end
%% CONVERT DATA POINTS
imageFRC = zeros(binCount, 2, zSlices);
binThresholds = zeros(binCount, 3, zSlices);
for z = 1:zSlices
    for b = 1:binCount
        if binCounter(b, 1, z) > 0 && binValues(b, 1, z) > 0 && binValues(b, 2, z) > 0
            imageFRC(b, 2, z) = binValues(b, 3, z) / sqrt(binValues(b, 1, z) * binValues(b, 2, z));
        else
            imageFRC(b, 2, z) = NaN;
        end
        imageFRC(b, 1, z) = ((b * maxRadius) / binCount) / (max(height, width) * pixelSize);
        if binCounter(b, 1, z) > 0
            binThresholds(b, 1, z) = 3 / sqrt(binCounter(b, 1, z) / 2);
            binThresholds(b, 2, z) = (0.2071 + 1.9102 / sqrt(binCounter(b, 1, z))) / (1.2071 + 0.9102 / sqrt(binCounter(b, 1, z)));
        else
            binThresholds(b, :, z) = NaN;
        end
    end
end
%% CALCULATE RESOLUTIONS
resolutions = zeros(3, 2, zSlices);
resolutionSum = zeros(3);
imageFRCSum = zeros(binCount, 2);
for z = 1:zSlices
    imageFRC(:, 2, z) = smoothdata(imageFRC(:, 2, z), 'movmean', max(round(binCount * binPercent), window));
    binThresholds(:, 3, z) = 1/7;
    fprintf("Z-Slice = %d\n", z);
    for t = 1:3
        crossingIndex = find(diff(sign(imageFRC(:, 2, z) - binThresholds(:, t, z))) < 0, 1, 'first');
        if isempty(crossingIndex) || isnan(crossingIndex)
            resolutions(t, :, z) = [NaN, NaN];
        else
            aboveSF = imageFRC(crossingIndex, 1, z);
            belowSF = imageFRC(crossingIndex + 1, 1, z);
            aboveFRC = imageFRC(crossingIndex, 2, z);
            belowFRC = imageFRC(crossingIndex + 1, 2, z);
            cutoffFrequency = aboveSF + ((binThresholds(crossingIndex, t ,z) - aboveFRC) / (belowFRC - aboveFRC)) * (belowSF - aboveSF);
            resolutions(t, :, z) = [binThresholds(crossingIndex, t, z), 1000 / cutoffFrequency];
        end
    end
    fprintf("Thresholds:\n\t3-Sigma = %e\n", binThresholds(crossingIndex, 1, z));
    fprintf("\tHalf-Bit = %e\n", binThresholds(crossingIndex, 2, z));
    fprintf("\tFixed-1/7 = %e\n", binThresholds(crossingIndex, 3, z));
    fprintf("Resolutions:\n\t3-Sigma: %g nm\n", resolutions(1, 2, z));
    fprintf("\tHalf-Bit: %g nm\n", resolutions(2, 2, z));
    fprintf("\tFixed-1/7: %g nm\n", resolutions(3, 2, z));
    fprintf("==========================\n");
    resolutionSum(1) = resolutionSum(1) + resolutions(1, 2, z);
    resolutionSum(2) = resolutionSum(2) + resolutions(2, 2, z);
    resolutionSum(3) = resolutionSum(3) + resolutions(3, 2, z);
    imageFRCSum(:, 1) = imageFRCSum(:, 1) + imageFRC(:, 1, z);
    imageFRCSum(:, 2) = imageFRCSum(:, 2) + imageFRC(:, 2, z);
end
resolutionMean = zeros(3);
resolutionMean(1) = resolutionSum(1) / zSlices;
resolutionMean(2) = resolutionSum(2) / zSlices;
resolutionMean(3) = resolutionSum(3) / zSlices;
imageFRCMean = zeros(binCount, 2);
imageFRCMean(:, 1) = imageFRCSum(:, 1) / zSlices;
imageFRCMean(:, 2) = imageFRCSum(:, 2) / zSlices;
fprintf("Mean Resolutions:\n\t3-Sigma: %g nm\n", resolutionMean(1));
fprintf("\tHalf-Bit: %g nm\n", resolutionMean(2));
fprintf("\tFixed-1/7: %g nm\n", resolutionMean(3));
%% PLOT RESULTS
figure;
plot(imageFRCMean(:, 1), imageFRCMean(:, 2), 'k-', 'LineWidth', 1.5); 
hold on;
colors = lines(3);
cutoffMarkers = nan(3, 1);
thresholdLabels = ["3σ", "½-bit", "1/7"];
for t = 1:3
    sf = 1000 / resolutionMean(t);
    if ~isnan(sf) && ~isinf(sf) && sf > 0
        xline(sf, '--', sprintf('%s: %.1f nm', thresholdLabels(t), resolutionMean(t)), 'Color', colors(t,:), 'LabelHorizontalAlignment','left', 'LabelOrientation','horizontal', 'LineWidth', 1.5);
    end
end
xlabel('Spatial Frequency (μm^{-1})');
ylabel('Circularly-Averaged Fourier Ring Correlation');
title(sprintf('FRC Analysis'));
ylim([min(imageFRCMean(:, 2)) max(imageFRCMean(:, 2))]);
xlim([min(imageFRCMean(:, 1)) max(imageFRCMean(:, 1))]);
legend('FRC Curve', 'Location', 'northeastoutside');
grid on;
box on;
