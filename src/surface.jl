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
> read_radflux(input_file::String)

Output is a Dictionary with all default variables:
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

    ncvars = ARM.sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = ARM.retrieveVariables(input_file, ncvars, attrvars=attrib)

    return output
end
# ----/

# end of script.
