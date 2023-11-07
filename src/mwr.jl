# Part of ARMtools.jl
# Containing functions related to MWR products e.g. LOS, RET
#

# *************************************************
# Function to read MWR RET products
"""
Function to read data from MWR RET product from ARM.

```julia-repl
julia> mwr = getMWRData(mwr_file)
julia> mwr = getMWRData(mwr_file; addvars=["alt", "lon", "lat"])
julia> mwr = getMWRData(mwr_file; onlyvars=["tbsky23"], attrvars=["doi"])
```
WHERE:
* mwr\\_file::String full path to file name to read,
* addvars::Vector{String} list of additional NetCDF variables to read,
* onlyvars::Vector{String} read only the listed variables,
* attrvars::Vector{String} additionally read the listed attributes.

The default data fields are:
* :time
* :LWP
* :IWV

Alternative variables can be:
* lat => North latitude,
* lon => East longitude,
* alt => altitude above mean sea level (instrument level),
* surface\\_vapor\\_pres =>  Intrument level vapor pressure,
* surface\\_pres =>  Instrument level pressure,
* surface\\_rh => Instrument level relative humidity,
* tbsky23 => Brightness Temperature at 23GHz [K],
* tbsky31 => Brightness Temperature at 31GHz [K],

Example of attributes:
* side_id
* doi

SUPPORTED PRODUCTS:
* MWR RET
* MWR LOS

"""
function getMWRData(in_file::String; addvars=[], onlyvars=[], attrvars=[], extras...)


    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    if isvariablein(in_file, "be_lwp")
        ncvars[:LWP] = "be_lwp";  # [g/m²]
        ncvars[:IWV] = "be_pwv";  # [cm]
        
        factor_lwp = 1f0;
        factor_iwv = 997f-2;  # [cm] -> [kg m⁻²]

        ncvars[:CBH] = "cloud_base_height"   # [km] AGL
        ncvars[:CLT] = "cloud_temp"   # [K]
        ncvars[:SFT] = "surface_temp"   # [K]
        
    elseif isvariablein(in_file, "liq")
        ncvars[:LWP] = "liq";  # [cm]
        ncvars[:IWV] = "vap";  # [cm]
        factor_lwp = 1f4;
        factor_iwv = 997f-2;  # [cm] -> [kg m⁻²]

        attvars = Dict(
            :elevation => "elevation",
            :azimuth => "azimuth",
            :wet => "wet_window"
        )
    else
        @warn "None know LWP variables found in $in_file"
        return nothing
    end
    
    # for MWR data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(in_file, ncvars, attrvars=attrvars)

    # converting units for LWP and IWV
    if haskey(output, :LWP)
        output[:LWP] *= factor_lwp
    end
    if haskey(output, :IWV)
        output[:IWV] *= factor_iwv
    end

    # MWR product RET has not info about elevation, azimuth nor wetness:
    if !haskey(output, :elevation)
        output[:elevation] = 90f0
    end
    if !haskey(output, :azimuth)
        output[:azimuth] = 0f0
    end
    if !haskey(output, :wet)
        output[:wet] = Int32(0)
    end
    
    return output
end
# ----/



# ---end of script.
