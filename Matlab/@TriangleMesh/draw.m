function h = draw(obj, varargin)
F=obj.F;
if size(obj.V,1)==2
    V=[obj.V;repmat(-1,1,length(obj.V))];
elseif size(obj.V,1)==3
    V=obj.V;
elseif isempty(obj.V)
    V=[];
end
if nargin==1
    h=patch;
    options={};
elseif nargin>=2
    if ishghandle(varargin{1})
        h=varargin{1};
        options=varargin(2:end);
    else
        h=patch;
        options=varargin(1:end);
    end
end

if isempty(options)
    options={'FaceColor',[0.7,0.7,0.7],'EdgeColor','black','CDataMapping','direct'};
end

set(h,'Vertices', V', 'Faces', F', options{:});
