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
    for (key_var, str_var) ∈ ncvars
        
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

        # *******************************
        # reading variable:
        tmp_var = ncin[str_var]
                
        # getting missing value attribute:
        fill_val = let vartype = eltype(skipmissing(tmp_var))
            tmp_fill = if haskey(tmp_var.attrib, "_FillValue")
                fillvalue(tmp_var)
            elseif haskey(tmp_var.attrib,  "missing_value")
                tmp_var.attrib["missing_value"]
            else
                -696969
            end
            vartype <: AbstractFloat ? vartype(NaN) : tmp_fill
        end
        
        tmp_var = let var_load=tmp_var[:]
            try
                if typeof(var_load)<:Array
                    NCDatasets.nomissing(var_load, fill_val)
                else
                    ismissing(var_load) ?  NaN32 : var_load
                end
            catch
                @warn "Variable $(str_var) seems empty or all ::Missing. Being skipped."
                continue
            end
        end
        
	# selecting only data within given limits (if provided in attributes):
        if haskey(ncin[str_var].attrib, "valid_min") && haskey(ncin[str_var].attrib, "valid_max")
            Vmin = ncin[str_var].attrib["valid_min"]
            Vmax = ncin[str_var].attrib["valid_max"]
            idx_out = (Vmin .≤ tmp_var .≤ Vmax) .|> !
            
            typeof(idx_out)<:BitArray ? (tmp_var[idx_out] .= NaN32) : (tmp_var = NaN32)
        end
        
        # filling the output variable:
        output[key_var] = tmp_var
    end
    
    # for global attributes:
    for (var, str_var) ∈ attrvars
        !haskey(ncin.attrib, str_var) && continue
        local tmp_out = ncin.attrib[str_var] #|> split

        
        if typeof(tmp_out) <: String
            tmp_out = map(split(tmp_out)) do x
                let tmp_type = any(contains.(x, (".", "-", "+"))) ? Float32 : Int32
                    q = tryparse(tmp_type, x)
                    isnothing(q) ? x : q
                end
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
# inquires instrument/product name by asking
# attribute datastream
"""
"""
function isinstrument(fname::String, name::String; prod_id=false)
    datastream = NCDataset(fname, "r") do ds
        haskey(ds.attrib, "datastream") ? ds.attrib["datastream"] : fname
    end
       
    return contains(datastream, lowercase(name))
end
# ----/

# ****************************************
# * request whether or not a varaible is in file
function isvariablein(fname::String, varname::String; attrib_flag=false)
    nc = NCDataset(fname, "r")
    isthere= attrib_flag ? haskey(nc.attrib, varname) : haskey(nc, varname)
    close(nc)
    return isthere
end

# ----/

# ******************************************
# * get data level from file name
function giveme_datalevel(fname::String)
    data_level = try
        nc = NCDataset(fname, "r")
        giveme_datalevel(nc)
    catch
        i1 = findfirst('.', fname)
        fname[i1+3] != '.' && error("$fname seems not to be ARM file?")
        fname[i1+1:i1+2]
    end
    if any(==(data_level), ["a1", "b1", "b2", "c0", "c1", "c2"])
        return data_level
    else
        @warn "ARM data level $(data_level) not recognized!"
        return nothing
    end
end
function giveme_datalevel(nc::NCDataset)
    
    haskey(nc.attrib, "data_level") && error("input dataset seems not to be ARM file?")
    return nc.attrib["data_level"]
end

# ****************************************
# * get file from pattern:
"""
Function to obtain list of files/single file which compraises with given criterium.

USAGE:
> files = getFilePattern(path::String, product::String, yy, mm ,dd; hh, submonth=false, fileext=".nc")

WHERE:
* path::String the path where to look files
* product::String Product name, then the files are scanned at "path/product"
* yy, mm, dd are the year, month and day, then the scanning is at "path/product/yy/*yymmdd*"
* hh (optional) is the hour, if specified the scanning follows "path/product/yy/*yymmdd.hh*"
* submonth (optional) is a flag to create sub-folder with month "path/product/yy/mm/*yymmdd*"
* fileext (optional) searchs for files with the extendion "path/product/yy/mm/*yymmdd*.fileext"

If it matchs multiple files the output will be a Vector{String} otherwise a String or nothing

"""
function getFilePattern(path::String, product::String, yy, mm, dd;
                        hh=nothing, submonth=false, fileext::String="")

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
    ofile = let tmp = occursin.(pattern, list_file) |> x->list_file[x]
        isempty(fileext) ? tmp : occursin.(fileext, tmp) |> x->tmp[x]
    end

    if typeof(ofile)<:Array && length(ofile)>1
        @warn "Multiple files match the pattern $(pattern) !!."
        return ofile
        
    elseif isempty(ofile)
        @warn "No files were found with the pattern $pattern !"
        return nothing
    else
        return ofile[1]
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
    return let ttx = extrema(ts) |> x-> x[1]:δts:x[2]
        [argmin(abs.(v.-ts)) for v in ttx]
    end
    #return extrema(ts) |> x-> x[1]:δts:x[2] |> x->findall(>(0), ts .∈ [x])
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

# *******************************************************************
# Surface related functions:
include("surface.jl")
# ----/

# *******************************************************************
# Infrared Temperature functions:
include("irt.jl")
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
