function [V,F,Fs] = offReader(filename)

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); Nv= str2num(a);
[a,str] = strtok(str); Nf= str2num(a);



[A,cnt] = fscanf(fid,'%f %f %f', 3*Nv);
if cnt~=3*Nv
    warning('Problem in reading vertices.');
end
V = reshape(A, 3, cnt/3);

% read Face 1  1088 480 1022
F=zeros(20,Nf);
Fs=zeros(1,Nf);
for i=1:Nf
    Fs(i)=fscanf(fid,'%d',1);
    F(1:Fs(i),i)=fscanf(fid,'%d',Fs(i))'+1;
end
F=F(1:max(Fs),:);
fclose(fid);



end

