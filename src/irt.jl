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
function getGNDIRTdata(in_file::Vector{String}; addvars=[], onlyvars=[], attrvars=[])
    dat_out = Dict{Symbol, Any}()
    catvar = Dict{Symbol, Union{Nothing, Int}}()
    ntime = -1

    getdim(x,n) = findall(==(n), size(x))
    
    foreach(in_file) do fn
        # reading single file:
        data = getGNDIRTdata(fn, addvars=addvars, onlyvars=onlyvars, attrvars=attrvars)

        if isempty(dat_out)
            dat_out = data
            catvar = let ntime = length(data[:time])
                tmp = [k=>getdim(v, ntime) for (k,v) ∈ data if typeof(v)<:Array]
                filter(p->!isnothing(p.second), tmp) |> Dict
            end
        else
            # merging every variable following time dimension v
            [dat_out[k] = cat(dat_out[k], data[k]; dims=v) for (k,v) ∈ catvar if !isempty(v)]
        end
    end
    return dat_out

end
# ----/

# end of file
