# Part of ARNtools.jl
# containing functions related to MWR products e.g. LOS, RET
#

# *************************************************
# Function to read MWR RET products
"""
Function getMWRData(mwr_file::String; addvars=[], onlyvars=[], attrvars=[])

This function read data from the Microwave Radiometer RET retrievals

The default data fields are:
* :time
* :LWP
* :IWV

Alternative variables can be:
* lat => North latitude
* lon => East longitude
* alt => altitude above mean sea level
* surface_vapor_pres => 
* surface_pres => 
* surface_rh =>
* tbsky23 => Brightness Temperature at 23GHz [K]
* tbsky31 => Brightness Temperature at 31GHz [K]

Attributes:
* side_id
* doi
"""
function getMWRData(mwr_file::String; addvars=[], onlyvars=[], attrvars=[])


    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    if isvariablein(mwr_file, "be_lwp")
        ncvars[:LWP] = "be_lwp";  # [g/m²]
        ncvars[:IWV] = "be_pwv";  # [cm]
        
        factor_lwp = 1f0;
        factor_iwv = 997f-2;  # [cm] -> [kg m⁻²]

        ncvars[:CBH] = "cloud_base_height"   # [km] AGL
        ncvars[:CLT] = "cloud_temp"   # [K]
        ncvars[:SFT] = "surface_temp"   # [K]
        
    elseif isvariablein(mwr_file, "liq")
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
        @error "None know LWP variables found in $mwr_file"
    end
    
    # for MWR data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(mwr_file, ncvars, attrvars=attrvars)

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
