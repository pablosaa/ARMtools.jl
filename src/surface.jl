# This file is part of the ARMtools set of function
# Here are the functions to read and analyze data related
# with surface properties measured by ARM


# ************************************************************
# Function to read SPECTRAL SURFACE ALBEDO
"""
Function get_SurfEnergyBudget(file_name::String)
"""
function read_SurfEnergyBudget(input_file::String; addvars=[], onlyvars=[], attrvars=[])
    # defaul netCDF variables to read from :
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    attrib = Dict()

    # For files type sebsE10.b1
    merge!(ncvars, Dict(:α => "albedo",
                        :surf_up_lw => "up_long",
                        :sky_dw_lw => "down_long",
                        :dw_sw => "down_short_hemisp",
                        :up_sw => "up_short_hemisp",
                        :net => "net_radiation",
                        :seb => "surface_energy_balance",
                        :Trad => "temp_net_radiometer",
                        :wet => "wetness",
                        )
           )

    merge!(attrib, Dict(:location => "facility_id",
                        :average_time => "averaging_interval",
                        )
           )

    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(input_file, ncvars, attrvars=attrib)
    return output
end

# **************************************************************
# Function to read RADFLUX data
"""
Function to read the data files for RADFLUX database
USAGE:
```julia-repl
julia> flx = read_radflux(input_file::String)
```
or, for a set of input files input as Vector{String}
```julia-repl
julia> flx = read_radflux(input_file::Vector{String})
julia> flx = read_radflux(base_path::String, dateofyear::Date)
```

Output is a Dictionary with standard varibles.

NOTE: The ARM rad flux data comes in daily files starting from 10:00UTC at the given
filename, and ending at 9:59UTC at the next day. Sometimes it is needed the data only
for a given Date, in that case use the input argument as:
```read_radflux("/data/nsa/", Date(year, month, day))```
which will output the data corresponding to the given day from 00:00UTC to 23:59UTC
found at the path "/data/nsa/yyyy/".
"""
function read_radflux(input_file::String; addvars=[], onlyvars=[], attrvars=[])
    # defaul netCDF variables to read from:
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    attrib = Dict()

    # For files type radfluxC1.b1
    merge!(ncvars, Dict(:up_lw => "upwelling_longwave",
                        :dw_lw => "downwelling_longwave",
                        :up_sw => "upwelling_shortwave",
                        :dw_sw => "downwelling_shortwave",
                        :sky_up_lw => "clearsky_upwelling_longwave",
                        :sky_dw_lw => "clearsky_downwelling_longwave",
                        :sky_up_sw => "clearsky_upwelling_shortwave",
                        :sky_dw_sw => "clearsky_downwelling_shortwave",
                        :transmisivity => "cloud_transmissivity_shortwave",
                        :cloud_TB => "cloud_radiating_temperature",
                        :sky_TB => "brightness_temperature",
                        :T_air => "air_temperature",
                        :cf_lw => "cloudfraction_longwave",
                        :cf_sw => "cloudfraction_shortwave",
                        :cos_za => "cosine_zenith",
                        )
           )

    merge!(attrib, Dict(:site => "site_id",
                        :location => "location_description",
                        :doi => "doi",
                        )
           )

    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(input_file, ncvars, attrvars=attrib)

    return output
end
function read_radflux(input_file::Vector{String}; addvars=[], onlyvars=[], attrvars=[])
    rad_out = Dict{Symbol, Any}()
    catvar = Dict{Symbol, Union{Nothing, Int}}()
    ntime = -1

    getdim(x,n) = findall(==(n), size(x))

    foreach(input_file) do fn
        # reading single file
        rad = read_radflux(fn, addvars=addvars, onlyvars=onlyvars, attrvars=attrvars)

        if isempty(rad_out)
            rad_out = rad
            catvar = let ntime=length(rad[:time])
                tmp = [k=>getdim(v, ntime) for (k,v) ∈ rad if typeof(v)<:Array]
                filter(p->!isnothing(p.second), tmp) |> Dict
            end
        else
            # merging every variable following time dimension v
            [rad_out[k] = cat(rad_out[k], rad[k]; dims=v) for (k,v) ∈ catvar if !isempty(v)]
        end
    end
    return rad_out
end
function read_radflux(path_name, thedate::Date; addvars=[], onlyvars=[], attrvars=[])
    DTs = [thedate - Day(1), thedate]
    thefiles = [ARMtools.getFilePattern(path_name,"", year(dt), month(dt), day(dt)) for dt ∈ DTs]
        
    rad_out = let rad=read_radflux(thefiles, addvars=addvars, onlyvars=onlyvars, attrvars=attrvars)
        idx = findall(==(thedate), Date.(rad[:time]))
        isempty(idx) && @warn("No data matched the $(thedate)")
        Dict(k=>typeof(v)<:Array && size(v,1)==length(rad[:time]) ? v[idx] : v for (k,v) ∈ rad)        
    end
    return rad_out
end
# ----/

# end of script.
