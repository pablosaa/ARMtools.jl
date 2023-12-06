# Part of ARMtools.jl
# Containing function related to Present Weather Disdrometer (PWD) sensors.
#
# ---
# TODO:
# [ ] include alternative variable list
# [ ] improve NaN replacement, right now invalid data is -9.223371f18
# [ ] include reading capability for Vector{Int} in retrievevariables()
# [ ] add option to include qc_ in the reading,
# [ ] add attribute information: e.g. platform_id, facility_id, etc.

# *****************************************************
# Function to read LD product
"""
This function returns data from Present Weather (PWD) sensors.

USAGE:
```julia-repl
julia> pwd= getPWDdata(pwd_file)
julia> pwd= getPWDdata(pwd_file; addvars=["pwd_mean_vis_10min"])
julia> pwd= getPWDdata(pwd_file; onlyvars=["PR", "rain"], attrvars=["site_id"])
```

WHERE:
* pwd\\_file::String full path of the PWD data file,
* addvars::Vector{String} list of netCDF variables to add to the default variables,

OUTPUT:
* pwd::Dict{Symbol, Vector} with the following default keys:
* :time
* :lat
* :lon
* :alt
* :vis   = "1 minute mean visibility"    # [m]
* :PR    = "1 minute mean precip_rate"   # [mm hr⁻¹]
* :Σrain = "cumulative rain"             # [mm]
* :Σsnow = "cumulative snow"             # [mm]
* :code  = "instantaneous weather code"  # [-]

Alternative variables can be:
*  => 
*  => 
*  => 

"""
function getPWDdata(in_file::String; addvars=[], onlyvars=[], attrvars=[])
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    if isvariablein(in_file, "pwd_precip_rate_mean_1min")
        #ncvars[:vis]  = "pwd_mean_vis_1min"           # [m]
        ncvars[:PR]   = "pwd_precip_rate_mean_1min"   # [mm hr⁻¹]
        ncvars[:qc_PR]   = "qc_pwd_precip_rate_mean_1min"   # [mm hr⁻¹]
        ncvars[:Σrain]= "pwd_cumul_rain"              # [mm]
        ncvars[:Σsnow]= "pwd_cumul_snow"              # [mm]
        #ncvars[:code] = "pwd_pw_code_inst"            # [-]
    else
        @warn "$in_file seems not to be a ARM PWD file"
        return nothing
    end

    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(in_file, ncvars, attrvars=attrvars)

    # NOTE: temporal solution for invalid data:
    foreach(output) do (k,v)
        if typeof(v)<:Vector && eltype(v)<:AbstractFloat
            ii = findall(<(0), v)
            output[k][ii] .= NaN
        end
    end
    
    return output
end
function getPWDdata(in_file::Vector{String}; addvars=[], onlyvars=[], attrvars=[])
    dat_out = Dict{Symbol, Any}()
    catvar = Dict{Symbol, Union{Nothing, Int}}()
    ntime = -1

    getdim(x,n) = findall(==(n), size(x))
    
    foreach(in_file) do fn
        # reading single file:
        data = getPWDdata(fn, addvars=addvars, onlyvars=onlyvars, attrvars=attrvars)

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
