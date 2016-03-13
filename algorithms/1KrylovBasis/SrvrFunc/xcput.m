function xcput(varargin)  
%eg:   CPU_time(time,   task,   repetitions Nr, order)
%--------------------------------------------------------------------------
narginchk(0,4);
SavingFlag = 0;
if nargin == 4
    t     = varargin{1}; 
    task  = varargin{2};    
    repeat= varargin{3};
    order = varargin{4};
elseif nargin == 3
    t     = varargin{1}; 
    task  = varargin{2};    
    repeat= varargin{3};
    order = 0;    
elseif nargin == 2
    t     = varargin{1};     
    task  = varargin{2};
    repeat= 1;
    order = 0;
elseif nargin == 1
    t     = varargin{1};  
    task  = 'other';
    repeat= 0;
    order = 0;
else
    SavingFlag = 1;
end
%--------------------------------------------------------------------------
persistent n CPUtim col_task col_t  col_repeat col_order header_of_table
if isempty(n);
   n=0;
    col_task  = 1;
    col_t     = 2;
    col_repeat= 3;
    col_order = 4;
    header_of_table ={'% #: ',  ...
                      'Task ',  ...
                      'Elapsed-Time',  ...
                      'For-Repitation  ',  ...
                      'Order-of-Sys.'};
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
if SavingFlag==1;
    dirFN = '!DO/zCPUtime.txt';
    DataSave(dirFN, header_of_table, CPUtim);
    clear CPUtim n;
    return;
end
%--------------------------------------------------------------------------
idx=[];
for j=1:n
    if and(strcmpi(task, CPUtim{j,col_task}), (order == CPUtim{j,col_order}))
       idx=[idx,j];
    end
end
% if n~=0
% [~,~,idx] = intersect(, cell2DataMat({}), 'rows')
% elselower()
%     idx=[];
% end

if isempty(idx);
    n=n+1;
    CPUtim{n, col_task } = lower(task);
    CPUtim{n, col_t } = t;
    CPUtim{n, col_repeat } =repeat; 
    CPUtim{n, col_order } =order;
elseif length(idx)==1;
    m=2;   CPUtim{idx,m} = CPUtim{idx,m} + t;
    m=m+1; CPUtim{idx,m} = CPUtim{idx,m} + repeat;  
else
    warning('Check it!');
end   

end  %func

  
function DataSave(filename, header ,DataMat)
fid=fopen(filename, 'at');
HL='%---------------------------------------------------------------------';
%--------------------------------------------------------------------------
% Date:
C=clock;
DD=strcat('Date: ',num2str(C(1)),'-', num2str(C(2)),'-', num2str(C(3)));
TT=strcat('Time: ',num2str(C(4)),':', num2str(C(5)),':', num2str(floor(C(4))));
fprintf(fid,'%% %s, \t%s\n',DD, TT);  

% fprintf(fid,'%% File-Name: %s\t', filename);
% DirFN = pwd;
% cd ..;
% tmp = cd(DirFN);
% fprintf(fid,'%% Simulation: %s', DirFN(length(tmp)+1:end) ); 
%--------------------------------------------------------------------------
       
       
%--------------------------------------------------------------------------  
% Header of Table: 
fprintf(fid,'%s', HL);
fprintf(fid,'\n%s\n');
for n=1:length(cellfun('length', header));
    fprintf(fid,'%-s\t',header{n});
end
fprintf(fid,'\n');
fprintf(fid,'%s\n', HL);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% DataMatrix:
n = length(cellfun('length', DataMat(:,1)));
totalTime = 0;
for i=1:n
    fprintf(fid, '%%%3d\t',i);
    m = length( DataMat(n,1:end) ); 
    %----------------------------------------------------
    fprintf(fid, '%s\t\t',DataMat{i,1});
    totalTime =  totalTime + DataMat{i,2};      
    for j=2:m
        fprintf(fid, '%10g\t',DataMat{i,j});
    end
    %----------------------------------------------------
        fprintf(fid, '\n');
end
fprintf(fid,'%s\n', HL);
fprintf(fid,'%% Total time:\t%g\n\n',  totalTime);
%--------------------------------------------------------------------------

fclose(fid);
fprintf('%% Total CPU-time:\t%g\n\n',  totalTime);
fclose all;
end %func    
    
    
    