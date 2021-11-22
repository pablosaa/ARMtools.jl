# This file is part of the ARMtools set of function
# Here are the functions to read and analyze data related
# with surface properties measured by ARM


# ************************************************************
# Function to read SPECTRAL SURFACE ALBEDO
"""
Function get_SurfEnergyBudget(file_name::String)
"""
function get_SurfEnergyBudget(input_file::String; addvars=[], onlyvars=[], attrvars=[])
    # defaul netCDF variables to read from KAZR:
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
# end of script.
