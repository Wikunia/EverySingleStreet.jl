# EverySingleStreet

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Wikunia.github.io/EverySingleStreet.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Wikunia.github.io/EverySingleStreet.jl/dev/)
[![Build Status](https://github.com/Wikunia/EverySingleStreet.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Wikunia/EverySingleStreet.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Wikunia/EverySingleStreet.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Wikunia/EverySingleStreet.jl)

This project is intended to be used to keep track of all the streets you've visited in your city. Do you maybe want to run every single street of your hometown? 

## About

This project provides tools for processing and analyzing street network data. It offers functionalities for tasks such as:

- Map matching: Mapping gps tracks to the most likely walked actual path based on the underlying street network
- Calculating statistics: 
    - How much of the city are already walked (taking into account only walkable roads)
    - How much of districts in the city are already walked (if district info is provided)

## Installation

This project is made using the [Julia programming language](https://julialang.org). To install this package run:

```julia
] add https://github.com/Wikunia/EverySingleStreet.jl

using EverySingleStreet
```
as it is currently not registered in the Julia registry.

## Rundown of the functionality

You can download the city of your choice like this:

```julia
path = "CLZ.json"
EverySingleStreet.download("Clausthal-Zellerfeld, Germany", path)
```

This will download the network of all streets and paths of ["Clausthal-Zellerfeld, Germany"](https://www.google.com/maps/place/38678+Clausthal-Zellerfeld,+Deutschland/@51.8064437,10.2692092,12z/data=!3m1!4b1!4m6!3m5!1s0x47a53bee117526c1:0x425ac6d94ac3f80!8m2!3d51.8080063!4d10.3407069!16zL20vMDE0cjMw?entry=ttu) and saves it as `CLZ.json`. 

Then we only want to keep the walkable paths of the network. This can be achieved by the following function:

```julia
EverySingleStreet.filter_walkable_json!(path)
```

We can create an internal representation for this map with the following command:

```julia
city_map = EverySingleStreet.parse_map(path);
```

This also computes some more information for faster processing later on. It works best for smaller cities or when run on a machine with a bit more RAM. Another variant is covered later for those smaller machines.

Now you can load in your gpx data and map it to your street network. The following provides you with some artificial way to create gpx data. You can read your available file however you want and put it in that format to make it work. 
Some options for standard [Strava](https://www.strava.com).

```julia
using Dates
using Geodesy
using TimeZones
path = [
    LLA(51.80665409406621, 10.335359255147063),
    LLA(51.806410, 10.335425),
    LLA(51.805925724425855, 10.335097770385982),
    LLA(51.80526071662234, 10.334263247682527),
    LLA(51.80494642844135, 10.334420156901205),
    LLA(51.8046968207022, 10.335259688289772),
    LLA(51.8043692602778, 10.335160446571944),
    LLA(51.80400106293405, 10.335266393838086)
]
times = [now()+i*Second(30) for i in 1:length(path)]
points = [EverySingleStreet.GPSPoint(p, ZonedDateTime(time, TimeZone("UTC"))) for (p, time) in zip(path, times)]

gpxfile = EverySingleStreet.GPXFile("test", points)
list_of_candidates = EverySingleStreet.map_path(city_map, gpxfile)
```

This maps the points to the walkable street network and returns a list of candidate points representing the most likely route. If there is a gap in between for which no possible path was found the list is split up into several parts.
Therefore it returns a vector of vectors of `Candidate`.

Let's have a look at the first candidate point:

```julia
first_candidate = candidates[1][1]
```

Each candidate point is in this format:

```julia
struct Candidate
    measured_point::GPSPoint
    lla::LLA
    way::Way
    way_is_reverse::Bool
    dist::Float64
    λ::Float64
end
```

Let's quickly have a look at each of those fields:

- measured_point: The actual recorded `GPSPoint` 
- lla: The latitude, longitude and altitude of the point on the street network. (Altitude doesn't contain reasonable values atm)
- way: A `Way` object that contains information about the street
- way_is_reverse: Whether the λ shows the distance rom the start or from the end of the way
- dist: Holds the information of the distance between the measured point and the point it was mapped to. `euclidean_distance(candidate.measured_point.pos, candidate.lla)`
- λ: Holds info of how far away this point is from the start (if `way_is_reverse` is `false`) or from the end (if `way_is_reverse` is `true`) along the way. This makes it simple to calculate which parts of each way are already visited

You can also get the information in a different format which might be more helpful in certain cases when tracking the achievement of how much one has walked already.

```julia
streetpaths = EverySingleStreet.calculate_streetpath("test", 1, candidates[1], city_map)
```

This returns a vector of `StreetPath` which consists of the following fields:

```julia
struct StreetPath
    name::String
    subpath_id::Int
    segments::Vector{StreetSegment}
end
```

The name and subpath_id are given in the previous function but the subpath_id can be increased if it was split up again into several parts. This can be true if there was no reasonable shortest path found between two candidate points as an example.

The most important field however is the `segments` field. Each segment is described in this struct:

```julia
struct StreetSegment
    from::Candidate
    to::Candidate
    function StreetSegment(from, to)
        @assert from.way.id == to.way.id
        @assert from.way_is_reverse == to.way_is_reverse
        new(from, to)
    end
end
```

So it simply consists of two `Candidate`s but for those it is asserted that they are both from the same street and in the same direction. 
For example if previously candidate 1 was in one street and candidate 2 in another this will create several segments such that this isn't the case anymore. For example by adding the end point of the street on which candidate 1 is on and the start point of the street where candidate 2 is on.

There is one more representation which is helpful when looking at several walks combined.

Let's have a look at a different way to load a city map which is useful for machines which can't hold a whole city network with precomputed fields in RAM. 

```julia
_, altona_map = EverySingleStreet.parse_no_graph_map("Altona.json");
```
Let's load a different map from a file for this (which is inside the `test/data` folder).

The `parse_no_graph_map` functionality parses the json file but doesn't create a network graph out of it and doesn't precompute shortest paths on the graph.

In general the idea of this approach is that whenever one wants to map a new walk onto that network extracts a local map of that walked area and computes the shortest paths on that instead of from the whole city.

We initialize a struct of walked parts of the city like this:
```julia
walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}());
```

and then load a walk into it which is stored in a json file as well. You can look at the format by checking this file which is part of the `test/data` folder as well.

```julia
altona_walk_path = joinpath("altona_walk.json");
```

This now matches the walk to the map 
```julia
mm_data = EverySingleStreet.map_matching(altona_walk_path, altona_map, walked_parts);
```

This returns a named tuple with the following fields:

```
added_kms
this_walked_road_km
walked_parts
```

It stores how many kilometers were added to the walked_parts, how many did you walked (on the road) and then an updated representation of all the walked parts of the city. 
In this case `added_kms` and `this_walked_road_km` will be the same as it was our first walk but later on `added` can be less than  `this_walked_road_km`.

Now let's have a look at the last part with yet another struct:

```julia
struct WalkedParts
    names::Dict{String, Vector{Int}} 
    ways::Dict{Int, WalkedWay}
end
```

The first struct matches street names to a list of integer ids and then the other matches those integer ids to walked ways. 

A walked way is described as following:

```julia
mutable struct WalkedWay
    way::Way
    parts::Vector{Tuple{Float64, Float64}}
end
```

A `Way` holds general information like the open street map id, all the nodes that describe the way. The name of the way, the total length and some more.

The `parts` section holds tuples of which parts of the way were walked by distance from the start point. For example 
```julia
(30.771703385452646, 67.50341022110413)
```
Means that one walked the part 30.7m to 67.5m from the start of the way. 

## Tracking districts

It's possible to get more fine grained statistics on district level as well if those are provided. One example of such a file you can again find in the `test/data` folder.

```julia
path = joinpath(@__DIR__, "..", "data", "Luebeck.json");
EverySingleStreet.download("Lübeck, Germany", path);
EverySingleStreet.filter_walkable_json!(path);
_, city_map = EverySingleStreet.parse_no_graph_map(path, joinpath(@__DIR__, "..", "data", "luebeck_districts.geojson"));
```

This gives a glimpse of those statistics:

```julia
walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}())
nt = EverySingleStreet.map_matching(joinpath(@__DIR__, "..", "data", "strava_luebeck.json"), city_map, walked_parts);
district_percentages = EverySingleStreet.get_walked_district_perc(city_map, collect(values(nt.walked_parts.ways)))
```

This will give the following result:

```julia
OrderedCollections.OrderedDict{Symbol, Float64} with 6 entries:
  :Innenstadt                 => 22.3847
  :None                       => 0.0
  Symbol("Sankt Jürgen")      => 0.0
  Symbol("Sankt Gertrud")     => 0.0
  Symbol("Sankt Lorenz Nord") => 0.0
  Symbol("Sankt Lorenz Süd")  => 0.0
```