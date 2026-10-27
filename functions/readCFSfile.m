function DataExtract = readCFSfile(file)
% Made by Fabien. I am not responsible if it dones not work for your
% application :)

% Input : file. Must be the path\file name of the file you want to convert
% from CFS to matlab
% Output : DataExtract. Structure that contains as many fields as there are
% channels in the CFS file you want to convert. The length of the data may
% not be the same in each field because each variable may not have the same
% frequency sampling.

READ=0;
DSVAR=1;
RL8=6;
[fhandle]=matcfs64c('cfsOpenFile',file,READ,0); % read only

[time,date,comment]=matcfs64c('cfsGetGenInfo',fhandle);
[channels,fileVars,DSVars,dataSections]=matcfs64c('cfsGetFileInfo',fhandle);

DataExtract.Acquisition_time = time ;
DataExtract.Acquisition_date = date ;
DataExtract.comment = comment ;
DataExtract.DataSections = dataSections ;

% FrameStart (seconds since Acquisition_time).
% Used later to match frames to the Brainsight neuronavigation samples.
% Left empty if not present.
DataExtract.FrameStart_s = [];
for v = 1:DSVars
    [~,varType,~,varDesc] = matcfs64c('cfsGetVarDesc',fhandle,v-1,DSVAR);
    if varType == RL8 && strcmp(strtrim(varDesc),'Start')
        frameStart = nan(dataSections,1);
        for i = 1:dataSections
            frameStart(i) = matcfs64c('cfsGetVarVal',fhandle,v-1,DSVAR,i,RL8);
        end
        DataExtract.FrameStart_s = frameStart;
        break
    end
end

if dataSections > 1
    dsVec=1:dataSections;
else
    dsVec=1;
end

for j=1:channels
    [channelName,yUnits,xUnits,dataType,dataKind,spacing,other]=matcfs64c('cfsGetFileChan',fhandle,j-1);
    channelName =  [ strrep(channelName,' ','') '_' num2str(j)] ;
    for i=1:length(dsVec)
       [startOffset,points,yScale,yOffset,xScale,xOffset]=matcfs64c('cfsGetDSChan',fhandle,j-1,dsVec(i));

        if startsWith(channelName,"ADC0") && i==1
            stim_offset=xOffset;                                            % to check potential offset with EMG
        end
        
        if i==1
            DataExtract.(channelName).dat = [] ;
            DataExtract.(channelName).FreqS = 1/xScale ;
            if startsWith(channelName,"EMG") || startsWith(channelName,"ADC0")
                DataExtract.(channelName).PointsPerFrame = zeros(dataSections,1);
            end
        end

        if points > 0
            startPt=0;
            [data]=matcfs64c('cfsGetChanData',fhandle,j-1,dsVec(i),startPt,points,dataType);
            data=(data*yScale)+yOffset;

            if startsWith(channelName,"EMG")                                % x offset correction
                nb_frm_offset=round((stim_offset-xOffset)/xScale);
                if nb_frm_offset<0                                          % If the EMG frame starts later than the stim, pad with zeros to realign
                    frm_offset=zeros(abs(nb_frm_offset),1);
                    data=[frm_offset ; data];
                elseif nb_frm_offset>0                                      % If the EMG frame starts earlier than the stim, trim to realign
                    data=data(nb_frm_offset+1:end);
                end
            end

            if startsWith(channelName,"EMG") || startsWith(channelName,"ADC0")
                % Number of samples this frame contributes to the concatenated
                % vector. Used later to map a sample index back to its frame of
                % origin, and to resample stim/EMG one section at a time.
                DataExtract.(channelName).PointsPerFrame(i) = numel(data);
            end

            DataExtract.(channelName).dat = [ DataExtract.(channelName).dat ; data ] ;
        end
    end
end
ret=matcfs64c('cfsCloseFile',fhandle); % close the file