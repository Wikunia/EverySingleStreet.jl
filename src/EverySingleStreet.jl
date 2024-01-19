module EverySingleStreet

using Accessors
using DataFrames
using DataStructures
using Dates
using Distributions
using FileIO
using Geodesy
using GeoJSON
import GeometryBasics: Point2
using Graphs
using HMMBase
using JLD2
using JSON3
using LightOSM
using LightXML: XMLDocument, save_file, create_root, new_child, add_text, set_attribute, set_attributes
using LinearAlgebra
using NearestNeighbors
using OrderedCollections
using PolygonInbounds
using Preferences
using ProgressMeter
using SparseArrays
using StaticGraphs
using Statistics
using TimeZones
using Unitful

const CANDIDATES_MAXIMUM_DISTANCE = @load_preference("CANDIDATES_MAXIMUM_DISTANCE", 50)
const CANDIDATES_FILTER_DISTANCE = @load_preference("CANDIDATES_FILTER_DISTANCE", 25)
const GPS_STD_DEV = @load_preference("GPS_STD_DEV", 15)
const LOCAL_MAP_PADDING = @load_preference("LOCAL_MAP_PADDING", 200)

include("gpx.jl")
include("types.jl")
include("utils.jl")
include("calculations.jl")
include("shortest_path.jl")
include("districts.jl")

end
