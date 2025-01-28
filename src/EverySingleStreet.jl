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
import LightXML
using LightXML: XMLDocument, save_file, create_root, new_child, add_text, set_attribute, set_attributes
using LinearAlgebra
import Luxor
import Luxor: @png, @svg
using NearestNeighbors
using OrderedCollections
using PolygonInbounds
using Preferences
using ProgressMeter
using SimpleWeightedGraphs
using SparseArrays
using StaticGraphs
using Statistics
using TimeZones
using Unitful

const DEFAULT_PREFERNECE_VALUE = Dict(
    "CANDIDATES_MAXIMUM_DISTANCE" => 50,
    "CANDIDATES_FILTER_DISTANCE" => 25,
    "GPS_STD_DEV" => 10,
    "LOCAL_MAP_PADDING" => 200,
    "EXTEND_WALKED_WAY_UP_TO" => 10,
    "MAX_FILL_CYCLE_LENGTH" => 100,
    "MIN_FILL_CYCLE_PERC" => 30,
    "MUL_NON_WALKABLE_ROAD" => 2.0,
    "MUL_WALKED_ROAD" => 4.0,
)

include("gpx.jl")
include("types.jl")
include("utils.jl")
include("calculations.jl")
include("shortest_path.jl")
include("districts.jl")
include("drawing.jl")
include("routing.jl")

end
