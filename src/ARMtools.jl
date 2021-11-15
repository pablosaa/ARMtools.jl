"""
Julia Package including functions to read and process data from the DOE ARM program.

The package is mainly focused to work with instruments located at the NSA facility in Utqiagvik, Alaska.
Nonetheless since most of the ARM instrumentation is standard, the package can also be used for other facilities.

(c) 2020, Pablo Saavedra Garfias
University of Leipzig
Faculty of Physics and Geosciences
LIM

See LICENSE
"""


module ARMtools

# *************************************************************
# Adding needed packages:
using NCDatasets
using Statistics
using Printf
using Wavelets
using Dates

# *************************************************************
## * Auxiliary functions:
## 1) Define variables to be read fron netCDF files
function sortVariables(defvars; onlyvars=[], addvars=[])

    outvars = Dict(:time=>"time")
    if ~isempty(onlyvars)
        # finding whether the variables are already as default:
        idx = map(x->in(x, onlyvars), values(defvars)) |> findall
        @assert ~isempty(idx) "selected only variables not present in pre-set"

        # converting only variables to Dict:
        thekeys = collect(keys(defvars))
        thevalues = collect(values(defvars))
        tmp = map(x->thekeys[x] => thevalues[x], idx) |> Dict
        #tmp = map(x->Symbol(uppercase(x))=>x, onlyvars[idx]) |> Dict
        outvars = merge(outvars, tmp)
    else
        outvars = defvars
    end
    
    if ~isempty(addvars)
        idx = map(x->findfirst(==(x), outvars), addvars) .|> isnothing |> findall
        tmp = map(x->Symbol(uppercase(x))=>x, addvars[idx]) |> Dict
        outvars = merge(outvars, tmp)
    end
    
    return outvars
end
# ----/

## 2) Read variables from netCDF file:
function retrieveVariables(ncfile::String, ncvars; attrvars=[])

    output = Dict()
    @assert isfile(ncfile) "reading $ncfile but not found!"
    ncin = NCDataset(ncfile)
    for var ∈ ncvars
        str_var = var[2]
        
        # checking for scaling attributes to change the variable values:
        scale_factor = missing
        if haskey(ncin[str_var].attrib, "scale_factor")
            scale_factor = ncin[str_var].attrib["scale_factor"]
        end

        add_offset = missing
        if haskey(ncin[str_var].attrib, "add_offset")
            add_offset = ncin[str_var].attrib["add_offset"]
        end

        # checking whether or not the variable is logarithmic:
        isvarlog = false
        if haskey(ncin[str_var].attrib, "units")
            unit_str = ncin[str_var].attrib["units"] |> lowercase
            isvarlog = map(x->contains(unit_str,x),
                           ("log", "logarithm", "log10", "dbz", "dbm","db")) |> any
        end

        # getting missing value attribute:
        tmp_var = ncin[str_var][:,:]
        
        if haskey(ncin[str_var].attrib, "missing_value")
            miss_val = ncin[str_var].attrib["missing_value"]
            type_var = eltype(tmp_var)
            tmp_var[tmp_var .≈ miss_val] .= type_var <: AbstractFloat ? NaN : fillvalue(type_var)

        end

        # in case the variable has attribute _FillValue instead:
        if haskey(ncin[str_var].attrib, "_FillValue")
            type_var = fillvalue(ncin[str_var])
            tmp_var = NCDatasets.nomissing(tmp_var, eltype(type_var) <: AbstractFloat ? NaN : type_var)
        end
        
        
        # filling the output variable:
        key_var = var[1]
        output[key_var] = tmp_var
    end
    
    # for global attributes:
    for (var, str_var) ∈ attrvars
        haskey(ncin.attrib, str_var) ? println(var) : continue
        tmp_var = ncin.attrib[str_var] |> split

        local tmp_out = map(tmp_var) do x
            let tmp_type = any(contains.(x, (".", "-", "+"))) ? Float32 : Int32
                q = tryparse(tmp_type, x)
                isnothing(q) ? x : q
            end
        end
          
        output[var] = length(tmp_out)>1 ? tmp_out : tmp_out[1]
    end
    close(ncin)
    return output
end
# ----/

# ****************************************
# HELPER FUNCTIONS
# ****************************************
#

# ****************************************
# * request whether or not a varaible is in file
function isvariablein(fname::String, varname::String)
    nc = NCDataset(fname, "r")
    isthere= haskey(nc, varname)
    close(nc)
    return isthere
end
# ----/

# ****************************************
# * get file from pattern:
"""
 Function getFilePattern(path::String, product::String, yy, mm ,dd; hh, submonth=false)

 retrieve file name pattern based on year, month, day to get
 a string to read.
"""
function getFilePattern(path::String, product::String, yy, mm, dd;
                        hh=nothing, submonth=false)

    yyyy_mm_dir = submonth ? @sprintf("%04d/%02d", yy, mm) : @sprintf("%04d", yy)
    base_dir = joinpath(path, product, yyyy_mm_dir)
        
    @assert isdir(base_dir) error("$base_dir seems it does not exist!")
    
    list_file = readdir(base_dir, join=true)
    if !isnothing(hh)
        pattern = @sprintf("%04d%02d%02d.%02d", yy, mm, dd, hh)
    else
        pattern = @sprintf("%04d%02d%02d", yy, mm, dd)
    end
    ofile = filter(x->all(occursin.(pattern, x)), list_file)

    if typeof(ofile)<:Array && length(ofile)>1
        @warn "Multiple files match the pattern $pattern, but the first one returned."
        return ofile[1]
        
    elseif isempty(ofile)
        @warn "No files were found with the pattern $pattern !"
        return nothing
    else
        return ofile
    end

end
# ----/

# ********************************************
# Return the indexes of subset for the time
# dimension given a desired time resolution
"""
Given a vector with elements type ts::Vector{DateTime}, this function
 returns the indexes that extract the elements of ts that match a given
 resolution, e.g. for instance the time variable from a netCDF file with 30
 seconds resolution, to be extracted at 1 minute (default), 30 seconds, or
 hourly, then use the function as following:

> idx_1min = index_at_time_resolution(ts)
> idx_30sec = index_at_time_resolution(ts, δts = Dates.Second(30))
> idx_1hour = index_at_time_resolution(ts, δts = Dates.Hour(1))

then
> ts[idx_1min]  # minute vector
> ts[idx_30sec] # 30 seconds vector
> ts[idx_1hour] # hourly vector

"""
function index_at_time_resolution(ts::Vector{DateTime}; δts=Dates.Minute(1))

    return extrema(ts) |> x-> x[1]:δts:x[2] |> x->findall(>(0), ts .∈ [x])
end
# ----/

# ********************************************
# ********************************************
# READING FUNCTIONS

# *******************************************************************
# Radio Sonde functions:
include("rawsonde.jl")

# ******************************************************************
# LIDARS functions:
include("lidars.jl")
# ----/

# *******************************************************************
# KAZR functions:
include("kazr.jl")
# ----/

# *******************************************************************
# MWR functions:
include("mwr.jl")
# ----/

# *******************************************************************
# NAV functions:
include("nav.jl")
# ----/

end # module


# Main file containing the package module
# See LICENSE

# ****************************************
# * Calculating netCDF files by factors
##function convert_factor_offset(var_in::AbstractFloat,
##                               scale_factor::AbstractFloat,
##                               add_offset::AbstractFloat,
##                               isvarlog::Bool)
##    var = var_in
##    if isvarlog
##        var = scale_factor*10f0^(var/10f0)
##        var += 10f0^(add_offset/10f0)
##        var = 10f0*log10(var)
##    else
##        var *= scale_factor
##        var += add_offset
##    end
##
##    return var
##end
##function convert_factor_offset(var::AbstractArray,
##                               scale_factor::eltype(var),
##                               add_offset::eltype(var),
##                               isvarlog::Bool)
##    vararray = convert_factor_offset.(var, scale_factor, add_offset, isvarlog)
##    return vararray
##end
# ----/
