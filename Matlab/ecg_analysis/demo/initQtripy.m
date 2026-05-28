function q = initQtripy(qtripy_path, qtriplot_exe_path)
    insert(py.sys.path, int32(0), qtripy_path);
    py.importlib.import_module('qtripy')
    
    q = py.qtripy.QTripy(qtriplot_exe_path); 
    q.begin(int32(1041));
    q.reset();
end