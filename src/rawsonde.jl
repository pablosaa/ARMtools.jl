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

# *******************************************************
"""
Function to convert DateTime vector into fraction of 24 hours
> T_hr = DateTime2HoursOfDay(DateTime(2004, 11, 30, 11, 30, 00) )
> T_hr = DateTime2HoursOfDay(T::Vector{DateTime})
"""
DateTime2HoursOfDay(T::DateTime) = hour(T) + minute(T)/60 + second(T)/3600
DateTime2HoursOfDay(T::Vector{DateTime}) = DateTime2HoursOfDay.(T)

#= *******************************************************
Auxiliary function to merge radiosonde with GND TIR data
for the first level as altitude zero by default
=#
"""
Function to merge surface Temperature to Radiosonde data.

USAGE:

julia> attach\\_Tₛ!(RS::Dict, Tₛ::Vector, time\\_Tₛ::Vector)

julia> attach\\_Tₛ!(RS::Dict, Tₛ::Vector, time\\_Tₛ::Vector; Hₛ=2)

julia> attach\\_Tₛ!(RS::Dict, IRT::Dict)

WHERE
* RS : the radiosonde data as output by getSondeData()
* Tₛ : Vector with the temperature in °C
* time_Tₛ : Vector with the time for Tₛ
* Hₛ : height in m above surface to attach Tₛ (default 0 m)
* IRT: the ground infrared thermomether data as output by getGNDURTdata()

OUTPUT:
A modified Radiosonde data Dictionary.
"""
function attach_Tₛ!(RS::Dict, Tₛ::Vector, time_Tₛ::Vector; Hₛ=0.0)

    Hₛ *= 1f-3   # converting m to km since RS[:height] is in km

    Hₛ ≥ RS[:height][1] && @warn "Hₛ surface is higher or equal to lowest RS level!"
    rs_thr = DateTime2HoursOfDay(RS[:time]);
    ir_thr = DateTime2HoursOfDay(time_Tₛ);

    nodes = (RS[:height], rs_thr);

    foreach(RS) do (k, V)
    
        T = typeof(V)
        X = T<:Vector ? T(undef, length(V)+1) : T(undef, (size(V).+(1,0))...)
        X[2:end, : ] = V
        
        if k == :time
            X = V
        elseif k == :T
            rs_extrap = LinearInterpolation(ir_thr, Tₛ, extrapolation_bc=Line())
            X[1,:] = rs_extrap(Hₛ, rs_thr)
        elseif k == :height
            X = vcat(Hₛ, V)
        else
            rs_extrap = LinearInterpolation(nodes, V, extrapolation_bc=Line());
            X[1,:] = rs_extrap.(Hₛ, rs_thr)
        
        end
        RS[k] = X
    end
end
function attach_Tₛ!(RS::Dict, IRT::Dict; Hₛ=0.0)
    !haskey(IRT, :time) && (@error "second argument IRT::Dict needs key :time")
    !haskey(IRT, :IRT)  && (@error "seocnd argument IRT::Dict needs key :IRT in K")
    Ts = IRT[:IRT] .- 273.15
    time_Ts = IRT[:time]
    ARMtools.attach_Tₛ!(RS, Ts, time_Ts; Hₛ=Hₛ)
end

# --- end of script.
