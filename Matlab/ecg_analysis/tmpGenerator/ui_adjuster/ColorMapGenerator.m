classdef ColorMapGenerator
    properties
        pygenerator
        signals_data
    end
    methods (Access = public)
        function obj = ColorMapGenerator()
            insert(py.sys.path, int32(0), 'C:\Users\Admin\Documents\Projects\ecg project\Scripts\Python\MapGenerator');
            py.importlib.import_module('map_generator')
    
            obj.pygenerator = py.map_generator.MapGenerator('C:\Users\Admin\Documents\Projects\ecg project\Scripts\Python\MapGenerator\Human_body.png');
            obj = obj.getBSM64points();
        end

        function img = BSM64Corr2Img(obj, corr)
            values = obj.corr2ImgCorr(corr);
            img = obj.generate_colormap(values);
        end

    end
    methods (Access = private)
        function map = generate_colormap(obj, values)
            xy_matrix = cell2mat(obj.signals_data(:, 2:3));
            py_points = py.list();
            for i = 1:size(xy_matrix, 1)
                py_points.append(py.tuple({int32(xy_matrix(i,1)), int32(xy_matrix(i,2))}));
            end
            
            py_values = py.list(values(:)');
            py_rgb_image = obj.pygenerator.generate_colormap(py_points, py_values, draw_points=true);
            map = uint8(py_rgb_image);
        end

        function obj = getBSM64points(obj)
            obj.signals_data = {
                % --- Kolumna 1 (X = 195) ---
                '10', 195, 186;
                '11', 195, 257;
                '12', 195, 329;
                '13', 195, 400;
            
                % --- Kolumna 2 (X = 213) ---
                '14', 213, 222;
                '15', 213, 293;
                '16', 213, 365;
            
                % --- Kolumna 3 (X = 231) ---
                '17', 231, 186;
                '18', 231, 222;
                % 'V1', 231, 257;
                '20', 231, 293;
                '21', 231, 329;
                '22', 231, 365;
                '23', 231, 400;
            
                % --- Kolumna 4 (X = 249) ---
                '24', 249, 186;
                '25', 249, 222;
                % 'V2', 249, 257;
                '27', 249, 293;
                '28', 249, 329;
                '29', 249, 365;
                '30', 249, 436;
            
                % --- Kolumna 5 (X = 266) ---
                '31', 266, 186;
                '32', 266, 222;
                % 'V3', 266, 257;
                '34', 266, 293;
                '35', 266, 329;
                '36', 266, 365;
                '37', 266, 400;
            
                % --- Kolumna 6 (X = 284) ---
                '38', 284, 186;
                '39', 284, 222;
                '40', 284, 257;
                % 'V4', 284, 293;
                '42', 284, 329;
                '43', 284, 365;
                '44', 284, 436;
            
                % --- Kolumna 7 (X = 302) ---
                '45', 302, 186;
                '46', 302, 222;
                '47', 302, 257;
                % 'V5', 302, 293;
                '49', 302, 329;
                '50', 302, 365;
                '51', 302, 400;
            
                % --- Kolumna 8 (X = 320) ---
                '52', 320, 222;
                '53', 320, 257;
                % 'V6', 320, 293;
                '55', 320, 365;
                '56', 320, 436;
            };
        end

        function values = corr2ImgCorr(obj, corr)
            all_ids_text = obj.signals_data(:, 1);
            numeric_ids = str2double(all_ids_text); % '10' -> 10,  'V1' -> NaN
            is_numeric_id = ~isnan(numeric_ids);
            clean_indices = numeric_ids(is_numeric_id);

            values = corr(clean_indices);
        end
    end
end