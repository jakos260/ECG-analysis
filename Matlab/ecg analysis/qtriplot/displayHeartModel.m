function displayHeartModel(ventri_ver, ventri_tri, thorax_ver, thorax_tri, time, focus)
    % outside the function run:
    % global qtriplot_port; qtriplot_port = 1041;
    
    qtriplot('reset')
    qmarker(focus, 'black', 5)
    % qtriplot('horizontal 2')
    % qtriplot('panel 1 1')
    % qtriplot(thorax_ver, thorax_tri)
    % qtriplot('trans 0.9')
    % qtriplot('edge y')
    qtriplot(ventri_ver, ventri_tri)
    qtriplot('trans 0.3')
    qtriplot(time)
    

    % qtriplot('panel 2 1')
    % qtriplot(ventri_ver, ventri_tri)
    % qtriplot('trans 0.3')
    % 
    % qtriplot('panel 2 1')
    % qtriplot(t_ref)
    % 
    % qtriplot('panel 2 1')
    % qmarker(focus, 'red', 5)

end