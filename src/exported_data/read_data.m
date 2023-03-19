function [reflection_data, transmission_data,source_data, t_data] = read_data(reflection_file, transmission_file, source_file, t_file)

    reflection_data = readmatrix(reflection_file);
    transmission_data = readmatrix(transmission_file);
    source_data = readmatrix(source_file);
    t_data = readmatrix(t_file);
    

end