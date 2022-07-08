function [threshold] = select_threshold_value(dataSource, windowSize)

sources = ['S2','L8','L7','S1','RADARSAT'];


switch dataSource
    case 'S2'
        windowSize_keys = {150,300,750, 900}; % NB: the 900m is downsampled resolution of 300m; just use same as 300m
        thresholds = [0.046 0.040 0.039, 0.040];
        threshold_container = containers.Map(windowSize_keys,thresholds);
    case 'L8'
        windowSize_keys = {150,300,750};
        thresholds = [0.048 0.051 0.065 ];
        threshold_container = containers.Map(windowSize_keys,thresholds);
    case 'L7' % NB: values are not representative: "undamaged" windows contain stripes
        windowSize_keys = {150,300,750,990};
%         thresholds = [0.045 0.058 0.071 0.073 ]; % based on 'undamaged' windows with heaavy striping
        thresholds = [0.027 0.032 0.037 0.035 ]; % based on 'undamaged' windows with less pronounced striping
        threshold_container = containers.Map(windowSize_keys,thresholds);
    case 'S1' 
        windowSize_keys = {150,300,750,1000};
        thresholds = [0.057 0.050 0.044 0.034];
        threshold_container = containers.Map(windowSize_keys,thresholds);
    case 'RADARSAT' 
        windowSize_keys = {1000};
%         thresholds = [0.047]; %[0.049]; % mean undamaged signals
        thresholds = [0.033]; % mean-std undamaged signals
        threshold_container = containers.Map(windowSize_keys,thresholds);    
end

if contains(sources, dataSource) & any(contains(string(windowSize_keys),string(windowSize)))
    threshold = threshold_container(windowSize);
else
    threshold = NaN;  
    disp(['[W] No threshold calculated for windowsize ' num2str(windowSize) ' and source ' dataSource ])
end


end