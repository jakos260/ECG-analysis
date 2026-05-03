function h =image3(A,fun)

if ndims(A)~=3
    error('alleen geschikt voor 3D images')
end
if nargin == 1
    fun = @imshow;
end

h.fig = figure;
setappdata(h.fig,'image',A)
setappdata(h.fig,'value',1)
h.image = feval(fun,A(:,:,1));

pos = get(gca,'position');
h.edit = uicontrol('style','edit','units','normalized','string','1');
set(h.edit,'position',[pos(1)+3/8*pos(3),pos(2)+pos(4)+.01, pos(3)/4, pos(4)/12],...
   'callback',{@adjust_slide,h});
set(h.image,'ButtonDownFcn',{@scroll,h})


function scroll(hObject,event,h)

A = getappdata(h.fig,'image');
val = getappdata(h.fig,'value');

if hObject == h.image
    click = get(h.fig,'SelectionType');
    switch click
        case {'normal','open'}
            val = val+1;
        otherwise
            val = val-1;
    end
end

if val>size(A,3)
    val = 1;
elseif val<1
    val = size(A,3);
end


set(h.image,'cdata',A(:,:,val))
set(h.edit,'string',num2str(val))
%set(h.title,'string',['Image number ' num2str(val)])
setappdata(h.fig,'value',val)

function adjust_slide(hObject,event,h)

str = get(h.edit,'string');
val = str2num(str);

if ~isempty(val)
    setappdata(h.fig,'value',val)
    scroll(hObject,event,h);
end
