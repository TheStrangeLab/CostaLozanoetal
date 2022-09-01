function [ eventClean, indFalseTriggers, TriggerValues_Steps, flagTriggerError ] = Get_Clean_Event_Structure( eventRaw )

% Clean event structure by detecting and deleting nonzero trigger values
% due to rising and falling edges of trigger signal

%% Get trigger values
TriggerValues = [eventRaw.value]';

%% Append zeros to the beginning and end of the trigger values array
% Consider multiple neighboring nonzero values also at the edges
TriggerValues = [0;TriggerValues;0];

%% Step 1: Locate doubles of nonzero triggers (2 nonzero triggers between 2 zeros)
% Get the location of zeros
indZeroLocations_1 = find(TriggerValues==0);
% Distance between zeros
indDistanceBetweenZeros_1 = indZeroLocations_1(2:end)-indZeroLocations_1(1:end-1);
% Zeros that have 2 nonzero values after them
indDoubleNonZerosRelative = indDistanceBetweenZeros_1==3;
indDoubleNonZeros = indZeroLocations_1(indDoubleNonZerosRelative);

%% Step 2 and 3 executed if there are doubles of nonzero triggers
cond_step_2_3 = ~isempty(indDoubleNonZeros);

if(cond_step_2_3) % if there are doubles of nonzero triggers
    
    %% Step 2: Detect false nonzero triggers on the falling edge
    % Neighboring nonzero triggers (Each line in the format: Zero - nonzero - nonzero - zero)
    DoubleNonZeroEntries_1 = TriggerValues([indDoubleNonZeros,indDoubleNonZeros+1,indDoubleNonZeros+2,indDoubleNonZeros+3]);
    if(length(indDoubleNonZeros)==1)
        DoubleNonZeroEntries_1 = DoubleNonZeroEntries_1';
    end
    % Nonzero triggers on the falling edge
    indDoubleNonZerosFallingRelative = DoubleNonZeroEntries_1(:,2)>=DoubleNonZeroEntries_1(:,3); % Relative to 'indDoubleNonZeros'
    indDoubleNonZerosFalling = indDoubleNonZeros(indDoubleNonZerosFallingRelative); % Relative to 'TriggerValues'
    % Store trigger values without nonzero triggers on the falling edge
    TriggerValues_2 = TriggerValues; % False nonzero triggers will be replaced by NaN
    TriggerValues_2(indDoubleNonZerosFalling+2) = NaN;
    
    %% Step 3: Detect false nonzero triggers on the rising edge
    % Nonzero triggers on the rising edge
    indDoubleNonZerosRisingRelative = DoubleNonZeroEntries_1(:,2)<DoubleNonZeroEntries_1(:,3); % Relative to 'indDoubleNonZeros'
    indDoubleNonZerosRising = indDoubleNonZeros(indDoubleNonZerosRisingRelative); % Relative to 'TriggerValues'
    % Store trigger values without nonzero triggers on the falling edge
    TriggerValues_3 = TriggerValues_2; % False nonzero triggers will be replaced by NaN
    TriggerValues_3(indDoubleNonZerosRising+1) = NaN;
    
    %%
else % if there are no doubles of nonzero triggers
    TriggerValues_2 = TriggerValues;
    TriggerValues_3 = TriggerValues;
end

%% Step 4: Locate triples of nonzero triggers (3 nonzero triggers between 2 zeros)
% Replace NaNs in 'TriggerValues_3' with zeros
TriggerValues_3_Zero = TriggerValues_3;
TriggerValues_3_Zero(isnan(TriggerValues_3_Zero)) = 0;
% Get the location of zeros
indZeroLocations_2 = find(TriggerValues_3_Zero==0);
% Distance between zeros
indDistanceBetweenZeros_2 = indZeroLocations_2(2:end)-indZeroLocations_2(1:end-1);
% Zeros that have 3 nonzero values after them
indTripleNonZerosRelative = indDistanceBetweenZeros_2==4;
indTripleNonZeros = indZeroLocations_2(indTripleNonZerosRelative);

%% Step 5: Keep the nonzero trigger in the center (the one with the largest value)
TripleNonZeroEntries_4 = TriggerValues_3([indTripleNonZeros,indTripleNonZeros+1,indTripleNonZeros+2,indTripleNonZeros+3,indTripleNonZeros+4]);
if(length(indTripleNonZeros)==1)
    TripleNonZeroEntries_4 = TripleNonZeroEntries_4';
end
% Non zero trigger in the center as the one with the largest value
indTripleNonZerosCenterMaxRelative = (TripleNonZeroEntries_4(:,2)<TripleNonZeroEntries_4(:,3))&(TripleNonZeroEntries_4(:,4)<=TripleNonZeroEntries_4(:,3)); % Relative to 'indTripleNonZeros'
indTripleNonZerosCenterMax = indTripleNonZeros(indTripleNonZerosCenterMaxRelative);
% Store trigger values without triples of nonzero values with the trigger in the center having the highest value
TriggerValues_4 = TriggerValues_3; % False nonzero triggers will be replaced by NaN
TriggerValues_4([indTripleNonZerosCenterMax+1;indTripleNonZerosCenterMax+3]) = NaN;

%% Step 6: Check again to locate multiple neighboring nonzero triggers
% Give an error if there are any
% Replace NaNs in 'TriggerValues_4' with zeros
TriggerValues_4_Zero = TriggerValues_4;
TriggerValues_4_Zero(isnan(TriggerValues_4_Zero)) = 0;
% Get the location of zeros
indZeroLocations_3 = find(TriggerValues_4_Zero==0);
% Distance between zeros
indDistanceBetweenZeros_3 = indZeroLocations_3(2:end)-indZeroLocations_3(1:end-1);
% Zeros that have more than one nonzero values after them
indMultipleNonZerosRelative = indDistanceBetweenZeros_3>2;
indMultipleNonZeros = indZeroLocations_3(indMultipleNonZerosRelative);

%% Give an error if there are multiple neighboring nonzero triggers
cond_error = ~isempty(find(indMultipleNonZerosRelative,1));
if(cond_error)
    flagTriggerError = 1;
else
    flagTriggerError = 0;
%     error('Error in trigger values')
end
% Consider other causes of false triggers if there is an error here
% Such as running out of battery

%% Delete values corresponding to zeros added to the beginning and the end of the trigger values array
TriggerValues(1) = [];
TriggerValues_2(1) = [];
TriggerValues_3(1) = [];
TriggerValues_4(1) = [];

TriggerValues(end) = [];
TriggerValues_2(end) = [];
TriggerValues_3(end) = [];
TriggerValues_4(end) = [];

%% Construct clean event structure
% False triggers are replaced by NaNs in 'TriggerValues_4'
indFalseTriggers = isnan(TriggerValues_4);
% Clean event structure
eventClean = eventRaw;
eventClean(indFalseTriggers) = [];

%% Initial trigger values and those after step 2, 3 and 5
TriggerValues_Steps = [TriggerValues,TriggerValues_2,TriggerValues_3,TriggerValues_4];


end