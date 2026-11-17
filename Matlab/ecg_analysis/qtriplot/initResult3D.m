function viewer = initResult3D(bins)
    if nargin < 1
        bins = 15;
    end

    eventPumpKey = 'Result3DEventPump';
    if isappdata(0, eventPumpKey)
        previousPump = getappdata(0, eventPumpKey);
        if isvalid(previousPump)
            stop(previousPump);
            delete(previousPump);
        end
        rmappdata(0, eventPumpKey);
    end

    insert(py.sys.path, int32(0), ...
        'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Python\qtripy');
    result3d = py.importlib.import_module('Result3D');

    viewer = result3d.Result3D(int32(bins));
    viewer.begin();
    viewer.reset();

    eventPump = timer( ...
        'ExecutionMode', 'fixedSpacing', ...
        'Period', 0.05, ...
        'TimerFcn', {@pumpResult3DEvents, viewer});
    setappdata(0, eventPumpKey, eventPump);
    start(eventPump);
end

function pumpResult3DEvents(eventPump, ~, viewer)
    if logical(viewer.process_events())
        return;
    end

    stop(eventPump);
    delete(eventPump);
    if isappdata(0, 'Result3DEventPump') && ...
            isequal(getappdata(0, 'Result3DEventPump'), eventPump)
        rmappdata(0, 'Result3DEventPump');
    end
end