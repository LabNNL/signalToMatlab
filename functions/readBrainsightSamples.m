function T = readBrainsightSamples(filepath)
%{
Parses the sample table of a Brainsight export (*.txt).
Returns a table with one row per sample:
Index, Target, Distance, TargetError, AngularError, TwistError and Time.
%}

lines=readlines(filepath);
lines=lines(strlength(lines)>0);

header_idx=find(startsWith(lines,"# Sample Name"),1);
if isempty(header_idx)
    error('No sample table found in : %s',filepath);
end

header=strtrim(split(lines(header_idx),sprintf('\t')));
header(1)=erase(header(1),"# ");

data_idx=header_idx+1;
while data_idx<=numel(lines) && ~startsWith(lines(data_idx),"#")
    data_idx=data_idx+1;
end
rows=split(lines(header_idx+1:data_idx-1),sprintf('\t'));

col=@(name) rows(:,header==name);

Index=double(col("Index"));
Target=col("Assoc. Target");
Distance=double(col("Dist. to Target"));
TargetError=double(col("Target Error"));
AngularError=double(col("Angular Error"));
TwistError=double(col("Twist Error"));
Time=timeofday(datetime(col("Time"),"InputFormat","HH:mm:ss.SSS"));

T=table(Index,Target,Distance,TargetError,AngularError,TwistError,Time);

end
