function viewer = initResult3D(bins)
    if nargin < 1
        bins = 15;
    end

    insert(py.sys.path, int32(0), ...
        'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Python\qtripy');
    result3d = py.importlib.import_module('Result3D');

    viewer = result3d.Result3D(int32(bins));
    viewer.begin();
    viewer.reset();
end