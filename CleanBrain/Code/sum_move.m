%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Function: This function estimates head motion based on output from FreeSurfer's fs_fast.

function [move] = sum_move(mov1, mov2, anfang, ende)

l = length(anfang);

if l > 1
    % Calculate motion metrics for multiple segments
    % (1) Between-slice motion or translation (in mm) 
    % (2) Within-slice up-down motion or translation (in mm) 
    % (3) Within-slice left-right motion or translation (in mm) 

    % Find the minimum and maximum values for the motion vectors in the specified ranges
    mini1 = min(mov1(anfang(1):ende(1), 1:3)); % Minimum for mov1 in specified range
    maxie1 = max(mov1(anfang(1):ende(1), 1:3)); % Maximum for mov1 in specified range

    mini2 = min(mov2(anfang(2):ende(2), 1:3)); % Minimum for mov2 in specified range
    maxie2 = max(mov2(anfang(2):ende(2), 1:3)); % Maximum for mov2 in specified range

    % Estimate the largest possible Euclidean distance between min and max motion vectors
    cart1 = maxie1 - mini1; % Motion span for mov1
    cart2 = maxie2 - mini2; % Motion span for mov2

    % Calculate the Euclidean norm (magnitude) of the motion spans
    cart1 = sqrt(sum(cart1.^2)); % Euclidean distance for mov1 in mm
    cart2 = sqrt(sum(cart2.^2)); % Euclidean distance for mov2 in mm

    % Extract the values corresponding to the volumes from fs_fast output
    movie1 = mov1(:, 4); % Volume data for mov1
    movie2 = mov2(:, 4); % Volume data for mov2

    % Define the segments for analysis based on onset and offset indices
    M1 = movie1(anfang(1):ende(1)); % Segment for mov1
    M2 = movie2(anfang(2):ende(2)); % Segment for mov2

    % Calculate the differences in volume data representing shifts in 3D space
    for i = 1:(length(M1) - 1)
        diff1(i) = M1(i + 1) - M1(i); % Difference for mov1
        diff2(i) = M2(i + 1) - M2(i); % Difference for mov2
    end

    % Aggregate motion metrics
    move.abs1 = sum(abs(diff1)); % Total absolute movement for mov1
    move.abs2 = sum(abs(diff2)); % Total absolute movement for mov2
    move.pow1 = sum((diff1).^2);  % Total squared movement for mov1
    move.pow2 = sum((diff2).^2);  % Total squared movement for mov2
    move.cart1 = cart1;           % Range of movement for mov1
    move.cart2 = cart2;           % Range of movement for mov2

else
    % Single segment processing when only one range is provided
    
    % Find minimum and maximum values for the motion vectors for mov1
    mini1 = min(mov1(anfang(1):ende(1), 1:3)); % Minimum for mov1 in specified range
    maxie1 = max(mov1(anfang(1):ende(1), 1:3)); % Maximum for mov1 in specified range
    cart1 = maxie1 - mini1; % Motion span for mov1

    % Extract the volume data for mov1
    movie1 = mov1(:, 4); % Volume data for mov1
    M1 = movie1(anfang(1):ende(1)); % Segment for mov1

    % Calculate the differences in volume data for the single segment
    for i = 1:(length(M1) - 1)
        diff1(i) = M1(i + 1) - M1(i); % Difference for mov1
    end

    % Aggregate motion metrics
    move.abs1 = sum(abs(diff1)); % Total absolute movement for mov1
    move.abs2 = nan;             % No second segment, set to NaN
    move.pow1 = sum((diff1).^2); % Total squared movement for mov1
    move.pow2 = nan;             % No second segment, set to NaN
    move.cart1 = cart1;          % Range of movement for mov1
    move.cart2 = nan;            % No second segment, set to NaN
    
end
