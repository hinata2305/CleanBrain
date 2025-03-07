% This code may not be used for medical applications. 
% ©Koten Schüppen 

% This code extracts reaction time data, identifies relevant time course beginnings and endings.

function [testA,testB,encode1A,encode2A,encode1B,encode2B,reactST,reactSR,reactMT,reactMR,respST,respSR,respMT,respMR] = extract_rt_rel_response_filter(subject)
  
    % Construct subject-specific log file name for subject 'a'
    subjectnameA = [subject 'a'];

    % Import log file for subject 'a'
    [out1, out2] = importPresentationLog([subjectnameA '-fage.log']);

    % Extract relevant fields from the imported log
    res = {out1.event_type};  % Event types from the log
    state = {out1.code};      % Corresponding state codes
    duration = {out1.ttime};  % Duration of each event

    % Identify indices for Stroop, Sternberg spatial and verbal memory tasks
    stroop = find(strcmp('Stroop Sternberg', state));
    mvis = find(strcmp('Retrieval Sternberg Symbols', state));
    mverb = find(strcmp('Retrieval Sternberg Letters', state));
    
    % Combine all responses that are waiting for user interaction
    wait_response = [stroop, mvis, mverb];

    % Initialize response and reaction time arrays
    response = nan(length(wait_response), 1);  % Preallocate for efficiency
    rtime = nan(length(wait_response), 1);

    % Loop over the indices of responses waiting from the log
    for ii = 1:length(wait_response)
        % Search for 'Response' event following the current wait_response
        tempr = find(strcmp(res(wait_response(ii) + 1 : wait_response(ii) + 6), 'Response'));

        if isempty(tempr)
            % If no response found, assign NaN
            response(ii) = nan;
            rtime(ii) = nan;
        elseif tempr < 4
            % If response is found in a valid index
            response(ii) = str2num(state{wait_response(ii) + tempr(1)});
            rtime(ii) = duration{wait_response(ii) + tempr(1)};
        else
            % If more than three indices away, adjust reaction time
            response(ii) = str2num(state{wait_response(ii) + tempr(1)});
            rtime(ii) = (duration{wait_response(ii) + tempr(1)}) + 25000; % Adding 25 seconds for delay
        end
    end

    % Extract reaction times and responses for spatial and verbal tasks for subject 'a'
    reactST = rtime(1:48);  % Reaction times for spatial task
    reactMT = rtime(49:96); % Reaction times for verbal task
    respST = response(1:48);  % Responses for spatial task
    respMT = response(49:96);  % Responses for verbal task

    % Identify the start and end indices for encoding time periods
    z1 = find(strcmp('Encoding Sternberg Symbols', state));
    z2 = find(strcmp('Encoding Sternberg Letters', state));
          
    e1 = find(strcmp('Instruction Sternberg Letters', state));
    e2 = find(strcmp('Instruction Stroop', state));
        
    % Find fixation cross event after retrieval
    r = find(strcmp('Fixation Cross after Retrieval', state));

    % Calculate the timing reference point from the first event in the log
    firstpulse = out1(1).time;

    % Calculate the encoding periods relative to the first pulse
    encode1A = round(([out1(z1).time] - firstpulse) ./ 12400);
    encode2A = round(([out1(z2).time] - firstpulse) ./ 12400);
      
    % Calculate end timings for the defined periods
    end1 = round(([out1(e1).time] - firstpulse) ./ 12400);
    end2 = round(([out1(e2).time] - firstpulse) ./ 12400);
    
    % Validate encoding timings based on expected values
    testA = (end1 - encode1A(1) == 487) & (end2 - encode2A(1) == 490);
       
    % Clear temporary variables
    clear rind res state duration

    % Repeat the process for subject 'b'
    subjectnameB = [subject 'b'];
    [out1, out2] = importPresentationLog([subjectnameB '-fage.log']);
  
    res = {out1.event_type};
    state = {out1.code};
    duration = {out1.ttime};
    
    % Identify indices for Stroop, Sternberg spatial and verbal memory tasks for subject 'b'
    stroop = find(strcmp('Stroop Sternberg', state));
    mvis = find(strcmp('Retrieval Sternberg Symbols', state));
    mverb = find(strcmp('Retrieval Sternberg Letters', state));
    
    wait_response = [stroop, mvis, mverb];

    % Initialize response and reaction time arrays for subject 'b'
    response = nan(length(wait_response), 1); 
    rtime = nan(length(wait_response), 1);

    % Loop over the indices of responses waiting from the log for subject 'b'
    for ii = 1:length(wait_response)
        % Search for 'Response' event
        tempr = find(strcmp(res(wait_response(ii) + 1 : wait_response(ii) + 6), 'Response'));

        if isempty(tempr)
            response(ii) = nan;  % No response found
            rtime(ii) = nan;     % No reaction time recorded
        elseif tempr < 4
            response(ii) = str2num(state{wait_response(ii) + tempr(1)});
            rtime(ii) = duration{wait_response(ii) + tempr(1)};
        else
            response(ii) = str2num(state{wait_response(ii) + tempr(1)});
            rtime(ii) = (duration{wait_response(ii) + tempr(1)}) + 25000; % Adjust for delay
        end
    end

    % Extract reaction times and responses for subject 'b'
    reactSR = rtime(1:48);  % Reaction times for spatial task
    reactMR = rtime(49:96); % Reaction times for verbal task
    respSR = response(1:48);  % Responses for spatial task
    respMR = response(49:96);  % Responses for verbal task

    % Identify the start and end indices for encoding time periods for subject 'b'
    z1 = find(strcmp('Encoding Sternberg Symbols', state));
    z2 = find(strcmp('Encoding Sternberg Letters', state));
      
    e1 = find(strcmp('Instruction Sternberg Letters', state));
    e2 = find(strcmp('Instruction Stroop', state));
        
    % Calculate the timing reference point for subject 'b'
    firstpulse = out1(1).time;

    % Calculate the encoding periods for subject 'b'
    encode1B = round(([out1(z1).time] - firstpulse) ./ 12400);
    encode2B = round(([out1(z2).time] - firstpulse) ./ 12400);
    
    % Calculate end timings for subject 'b'
    end1 = round(([out1(e1).time] - firstpulse) ./ 12400);
    end2 = round(([out1(e2).time] - firstpulse) ./ 12400);
    
    % Validate encoding timings for subject 'b'
    testB = (end1 - encode1B(1) == 487) & (end2 - encode2B(1) == 490);
   
end
