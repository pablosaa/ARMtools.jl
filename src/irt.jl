# Part of ARMtools.jl
# Containing function related to Infrared Radiometer Temperature
#

# *****************************************************
# Function to read GNDIRT product

"""
Function getGNDIRTdata(irt_file::String; addvars=[], onlyvars=[], attrvars=[])

This function returns data from Infrared Surface Temperature.

The default data fields are:
* :time
* :IRT
* :σIRT

Alternative variables can be:
* lat
* lon
* alt
* ref_ir_temp => reference internal IR temperature
* sfc_ir_temp_max => maximum IRT
* sfc_ir_temp_min => minimum IRT

"""
function getGNDIRTdata(in_file::String; addvars=[], onlyvars=[], attrvars=[])
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    if isvariablein(in_file, "sfc_ir_temp")
        ncvars[:IRT] = "sfc_ir_temp"
        ncvars[:σIRT] = "sfc_ir_temp_std"
    else
        @warn "$irt_file seems not to be a IRT file"
        return nothing
    end

    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(in_file, ncvars, attrvars=attrvars)

    return output
end
# ----/

# end of file
