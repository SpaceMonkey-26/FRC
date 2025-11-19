clc; clear;
%% INPUTS
inputImageLocation1 = "path1";
inputImageLocation2 = "path2";
pixelSize = 0.065; % Pixel size: micrometers per pixel
percentHighSF = 0.35; % Percentage of high spatial frequencies to estimate noise: 0.1-0.35
binPercent = 0.03; % Percentage of total bins for smoothing: 0.02-0.05
window = 7; % Minimum window size for smoothing: 3, 5, or 7
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
            bin = min((max(floor((radius / maxRadius) * binCount), 1) - 1), binCount - 1);
            binValues(bin + 1, 1, z) = binValues(bin + 1, 1, z) + imagePSD1(y, x, z);
            binValues(bin + 1, 2, z) = binValues(bin + 1, 2, z) + imagePSD2(y, x, z);
            binValues(bin + 1, 3, z) = binValues(bin + 1, 3, z) + imageCC(y, x, z);
            binCounter(bin + 1, 1, z) = binCounter(bin + 1, 1, z) + 1;
        end
    end
end
%% CONVERT DATA POINTS
imageFRC = zeros(binCount, 2, zSlices);
sortedFRC = zeros(height, width, zSlices);
binThresholds = zeros(binCount, 2, zSlices);
for z = 1:zSlices
    for b = 1:binCount
        if all(binValues(b, :, z) > 0)
            imageFRC(b, 2, z) = binValues(b, 3, z) / sqrt(binValues(b, 1, z) * binValues(b, 2, z));
        else
            imageFRC(b, 2, z) = NaN;
        end
        imageFRC(b, 1, z) = ((b * maxRadius) / binCount) / (max(height, width) * pixelSize);
        if binCounter(b, 1, z) > 0
            % binThresholds(b, 1, z) = (((sqrt(2) - 1) / 2) + (((sqrt(2) - 1) * sqrt(2)) / sqrt(binCounter(b, 1, z)))) / (1 + ((sqrt(2) - 1) / 2) + (((sqrt(2) - 1) * sqrt(2)) / sqrt(binCounter(b, 1, z))));
            % binThresholds(b, 2, z) = (((1 / 6) / 2) + (((1 / 6) * sqrt(2)) / sqrt(binCounter(b, 1, z)))) / (1 + ((1 / 6) / 2) + (((1 / 6) * sqrt(2)) / sqrt(binCounter(b, 1, z))));
        else
            binThresholds(b, :, z) = NaN;
        end
    end
    sortedFRC(:, :, z) = sortrows(imageFRC(:, :, z), -2);
end
%% CALCULATE RESOLUTIONS
resolutions = zeros(5, 2, zSlices);
thresholds = zeros(5);
for z = 1:zSlices
    imageFRC(:, 2, z) = smoothdata(imageFRC(:, 2, z), 'movmean', max(round(binCount * binPercent), window));
    sortedFRC(:, 2, z) = smoothdata(sortedFRC(:, 2, z), 'movmean', max(round(binCount * binPercent), window));
    noiseRegion = sortedFRC((floor(binCount * (1 - percentHighSF)) + 1):end, 2, z);
    mu = mean(noiseRegion);
    sigma = std(detrend(noiseRegion)); 
    thresholds(1) = mu + (3 * sigma);
    % thresholds(2) = ((1 - (1 / sqrt(2))) / (1 - (1 - (1 / sqrt(2))))) * (sigma^2);
    thresholds(3) = 1 - (1 / sqrt(2));
    % thresholds(4) = ((1 / 7) / (1 - (1 / 7))) * (sigma^2);
    thresholds(5) = 1 / 7;
    fprintf("Z-Slice = %d\n", z);
    fprintf("FRC Range: [%e, %e]\n", min(imageFRC(:, 2, z)), max(imageFRC(:, 2, z)));
    fprintf("Estimated Noise Level = %e\n", sigma);
    fprintf("Thresholds:\n\t3-Sigma = %e\n", thresholds(1));
    fprintf("\tHalf-Bit (Noise-dependent) = %e\n", thresholds(2));
    fprintf("\tHalf-Bit (Noise-independent) = %e\n", thresholds(3));
    fprintf("\tFixed-1/7 (Noise-dependent) = %e\n", thresholds(4));
    fprintf("\tFixed-1/7 (Noise-independent) = %e\n", thresholds(5));
    for t = 1:5
        % crossingIndex = find(diff(sign(imageFRC(:, 2, z) - thresholds(t))) < 0, 1, 'first');
        if isempty(crossingIndex)
            aboveFRC = find(imageFRC(:, 2, z) > thresholds(t));
            belowFRC = find(imageFRC(:, 2, z) < thresholds(t));
            if isempty(belowFRC) || isempty(aboveFRC)
                cutoffFrequency = NaN;
            else
                try
                    cutoffFrequency = interp1(imageFRC(:, 2, z), imageFRC(:, 1, z), thresholds(t), 'linear');
                catch
                    try
                        cutoffFrequency = interp1(imageFRC(:, 2, z), imageFRC(:, 1, z), thresholds(t), 'pchip');
                    catch
                        cutoffFrequency = NaN;
                    end
                end
            end
        else
            aboveSF = imageFRC(crossingIndex, 1, z);
            belowSF = imageFRC(crossingIndex + 1, 1, z);
            aboveFRC = imageFRC(crossingIndex, 2, z);
            belowFRC = imageFRC(crossingIndex + 1, 2, z);
            cutoffFrequency = aboveSF + ((thresholds(t) - aboveFRC) / (belowFRC - aboveFRC)) * (belowSF - aboveSF);
        end
        resolutions(t, :, z) = [thresholds(t), 1000 / cutoffFrequency];
    end
    fprintf("Resolutions:\n\t3-Sigma: %g nm\n", resolutions(1, 2, z));
    fprintf("\tHalf-Bit (Noise-dependent): %g nm\n", resolutions(2, 2, z));
    fprintf("\tHalf-Bit (Noise-independent): %g nm\n", resolutions(3, 2, z));
    fprintf("\tFixed-1/7 (Noise-dependent): %g nm\n", resolutions(4, 2, z));
    fprintf("\tFixed-1/7 (Noise-independent): %g nm\n", resolutions(5, 2, z));
    fprintf("=========\n");
%% PLOT RESULTS
    figure;
    plot(imageFRC(:, 1, z), imageFRC(:, 2, z), 'k-', 'LineWidth', 1.5); hold on;
    colors = lines(5);
    cutoffMarkers = nan(5, 1);
    thresholdLabels = ["3σ", "½-bit (ND)", "½-bit (NI)", "1/7 (ND)", "1/7 (NI)"];
    for t = 1:5
        yline(thresholds(t), '--', 'Color', colors(t, :), 'LineWidth', 1.2, 'Label', sprintf('%s', thresholdLabels(t), ' Threshold'), 'LabelHorizontalAlignment', 'left');
        sf = 1000 / resolutions(t, 2, z);
        if ~isnan(sf) && ~isinf(sf)
            cutoffMarkers(t) = sf;
            xline(sf, '-', 'Color', colors(t, :), 'LineStyle', '-', 'LineWidth', 1, 'Label', sprintf('%s: %.1f nm', thresholdLabels(t), resolutions(t, 2, z)), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
        end
    end
    xlabel('Spatial Frequency (μm^{-1})');
    ylabel('Circularly-Averaged Fourier Ring Correlation');
    title(sprintf('FRC Analysis – Z-Slice %d', z));
    ylim([min(imageFRC(:, 2, z)) max(imageFRC(:, 2, z))]);
    xlim([min(imageFRC(:, 1, z)) max(imageFRC(:, 1, z))]);
    legend('FRC Curve', 'Location', 'northeastoutside');
    grid on;
    box on;
end
