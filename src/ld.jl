# Part of ARMtools.jl
# Containing function related to Light sensor: e.g. Disdrometers
#
# ---
# TODO:
# [ ] include alternative variable list
# [ ] convert LWC to [g m⁻³]
# [ ] add attribute information: e.g. platform_id, facility_id, etc.


# *****************************************************
# Function to read LD product

"""
Function getLDdata(ld_file::String; addvars=[], onlyvars=[], attrvars=[])

This function returns data from Ligth sensors e.g. Disdrometers

The default data fields are:
* :time
* :lat
* :lon
* :alt
* :LWC = "liquid_water_content"    # [mm³ m⁻³]
* :PR = "precip_rate"              # [mm hr⁻¹]
* :D₀ = "median_volume_diameter"   # [mm]
* :Nₙ = "number_density_drops"     # [m⁻³ mm⁻¹]
* :Dᵢ = "particle_size"            # [mm]
* :Vₜ = "fall_velocity_calculated" # [m s⁻¹]

Alternative variables can be:
*  => 
*  => 
*  => 

"""
function getLDdata(in_file::String; addvars=[], onlyvars=[], attrvars=[])
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    if isvariablein(in_file, "raw_fall_velocity")
        ncvars[:LWC] = "liquid_water_content"    # [mm³ m⁻³]
        ncvars[:PR] = "precip_rate"              # [mm hr⁻¹]
        ncvars[:D₀] = "median_volume_diameter"   # [mm]
        ncvars[:Nₙ] = "number_density_drops"     # [m⁻³ mm⁻¹]
        ncvars[:Dᵢ] = "particle_size"            # [mm]
        ncvars[:Vₜ] = "fall_velocity_calculated" # [m s⁻¹]
    else
        @warn "$in_file seems not to be a ARM LD file"
        return nothing
    end

    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(in_file, ncvars, attrvars=attrvars)

    return output
end
function getLDdata(in_file::Vector{String}; addvars=[], onlyvars=[], attrvars=[])
    dat_out = Dict{Symbol, Any}()
    catvar = Dict{Symbol, Union{Nothing, Int}}()
    ntime = -1

    getdim(x,n) = findall(==(n), size(x))
    
    foreach(in_file) do fn
        # reading single file:
        data = getLDdata(fn, addvars=addvars, onlyvars=onlyvars, attrvars=attrvars)

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
