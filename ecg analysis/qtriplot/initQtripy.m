function q = initQtripy()
    insert(py.sys.path, int32(0), 'C:\Users\Admin\Documents\ecg project\qtripy');
    py.importlib.import_module('qtripy')
    
    q = py.qtripy.QTripy('C:\Program Files\qtriplot\qtriplot.exe');
    q.begin(int32(1041));
    q.reset();
end