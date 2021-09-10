# Part of ARMtools.jl
# Function relted to radio sondes

# ********************************************
# * Read INTERPOLATE radiosonde tools:
"""
Function getSondeData(file_name::String)

This will read the ARM datafile 'file_name.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :height
* :T
* :Pa
* :RH
* :qv
* :U
* :V
* :WSPD
* :WDIR
* :θ

"""
function getSondeData(sonde_file::String; addvars=[], onlyvars=[] )

    # defaul netCDF variables to read from Radiosonde:
    ncvars = Dict(:time=>"time",
                  :height=>"height",
                  :T=>"temp",
                  :Pa=>"bar_pres",
                  :RH=>"rh",
                  :qv=>"sh",
                  :U=>"u_wind",
                  :V=>"v_wind",
                  :WSPD=>"wspd",
                  :WDIR=>"wdir",
                  :θ=>"potential_temp")
    # for INTERPOLATE RS data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars)

    return output
end
# ----/

# --- end of script.
