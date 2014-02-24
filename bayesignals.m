% bayesignals.m
%
function waveform = bayesignals(DataStruct)

switch DataStruct.name
    case 'ringdown'
        % --- Time domain ringdown
            time=DataStruct.time;

            nonzeroidx=DataStruct.time>DataStruct.t0;
            waveform = zeros(1,length(time));
            waveform(nonzeroidx) = DataStruct.amp * ...
                cos(2*pi*DataStruct.freq*(DataStruct.time(nonzeroidx) - DataStruct.t0)) .* ...
                exp(-(DataStruct.time(nonzeroidx)-DataStruct.t0)/DataStruct.tau);
           
    otherwise
        disp('unknown signal')
        
end