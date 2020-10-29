%function to read the constants file in JSON and adapted to MAT
%inputs
function cst_p=read_CST_json(cst_file)
    
    struct = loadjson(cst_file);
    
    cst_p.semi_minor_axis_cst = struct.semi_minor_axis_cst.value;   
    cst_p.sec_in_day_cst = struct.sec_in_day_cst.value;
    cst_p.c_cst = struct.c_cst.value;
    cst_p.pi_cst = struct.pi_cst.value;
    cst_p.flat_coeff_cst = struct.flat_coeff_cst.value;
    cst_p.earth_radius_cst = struct.earth_radius_cst.value;
    cst_p.semi_major_axis_cst = struct.semi_major_axis_cst.value;

    clear struct;
end
