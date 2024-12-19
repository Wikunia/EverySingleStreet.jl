function add_district!(districts:: Vector{District}, name::String, geometry::GeoJSON.MultiPolygon)
    hpolygons = Vector{HolePolygon}()
    for polygon in geometry
        outer = polygon[1]
        holes = Vector{Vector{Tuple{Float32, Float32}}}()
        if length(polygon) > 1
            for hole in polygon[2:end]
                push!(holes, hole)
            end
        end
        push!(hpolygons, HolePolygon(outer, holes))
    end
    push!(districts, District(Symbol(name), hpolygons))
end

function add_district!(districts:: Vector{District}, name::String, geometry::GeoJSON.Polygon)
    hpolygons = Vector{HolePolygon}()
    for polygon in geometry.coordinates
        holes = Vector{Vector{Tuple{Float32, Float32}}}()
        push!(hpolygons, HolePolygon(polygon, holes))
    end
    push!(districts, District(Symbol(name), hpolygons))
end


"""
    get_districts(geojson_fpath)

Parse the districts in the given file which needs to have :geometry and :Stadtteil as properties.
Return a vector of [`District`](@ref)
"""
function get_districts(::Nothing)
    districts = Vector{District}()
    return districts
end

function get_districts(geojson_fpath)
    districts = GeoJSON.read(read(geojson_fpath))
    df_districts = DataFrame(districts) 
    districts = Vector{District}()
    for row in eachrow(df_districts)
        geometry = row[:geometry]
        add_district!(districts, row[:Stadtteil], geometry)
    end
    return districts
end

function points_in_polygon(points::Matrix, polygon)
    polygon_nodes = zeros(length(polygon), 2)
    for (i,node) in enumerate(polygon)
        # geojson have lon , lat order but we want lat lon order
        polygon_nodes[i,1] = node[2] # lat
        polygon_nodes[i,2] = node[1] # lon
    end
    inside_or_bound = inpoly2(points, polygon_nodes)
    inside = inside_or_bound[:, 1] .| inside_or_bound[:, 2]
    return Set(findall(inside))
end

function points_in_holepolygon(points::Matrix, holepolygon::HolePolygon)
    in_outer = points_in_polygon(points, holepolygon.outer)
    in_any_hole = Set{Int}()
    for hole in holepolygon.holes
        union!(in_any_hole, points_in_polygon(points, hole))
    end
    return setdiff(in_outer, in_any_hole)
end

function nodes_in_district(nodes::Vector{Node}, district::District)
    points = zeros(length(nodes), 2)
    for (i,node) in enumerate(nodes)
        points[i,1] = node.lat
        points[i,2] = node.lon
    end
    node_ids = Set{Int}()
    for holepolygon in district.polygons
        union!(node_ids, points_in_holepolygon(points, holepolygon))
    end
    return node_ids
end

function map_nodes_to_district(nodes::Vector{Node}, geojson_fpath::Nothing)
    node_district = Vector{Symbol}(undef, length(nodes))
    node_district .= :None
    return node_district
end

"""
    map_nodes_to_district(nodes::Vector{Node}, geojson_fpath)

Return a vector of district names (as Symbol) that points each [`Node`](@ref) to its respective district
"""
function map_nodes_to_district(nodes::Vector{Node}, geojson_fpath)
    districts = get_districts(geojson_fpath) 
    node_district = Vector{Symbol}(undef, length(nodes))
    node_district .= :None
    for district in districts 
        node_ids = collect(nodes_in_district(nodes, district))
        node_district[node_ids] .= district.name
    end
    return node_district
end

"""
    update_district_walked!(district_kms::AbstractDict{Symbol, Float64}, city_map::AbstractSimpleMap, walkedway::WalkedWay)

Update the total walked distance of each district after a [`WalkedWay`](@ref) is added to the map.

# Arguments
- `district_kms::AbstractDict{Symbol, Float64}`: A dictionary mapping district names (as Symbol) to their total walked distances in kilometers.
- `city_map::AbstractSimpleMap`: The city map object that stores nodes and edges.
- `walkedway::WalkedWay`: The walked way object added to the map.
"""
function update_district_walked!(district_kms::AbstractDict{Symbol, Float64}, city_map::AbstractSimpleMap, walkedway::WalkedWay, trans, rev_trans)
    way = walkedway.way
    nodes = way.nodes
    for part in walkedway.parts
        start_id = prev_idx(nodes, part[1])
        stop_id = prev_idx(nodes, part[2])+1
        # check if all ids are in the same district 
        district_names = Set{Symbol}()
        for i in start_id:stop_id
            node = nodes[i]
            haskey(city_map.osm_id_to_node_id, node.id) || continue
            internal_node_id = city_map.osm_id_to_node_id[node.id]
            district = city_map.nodes_to_district_name[internal_node_id]
            push!(district_names, district)
        end
        if length(district_names) == 1
            district = collect(district_names)[1]
            already_walked_in_district = get(district_kms, district, 0.0)
            dist = (part[2] - part[1]) / 1000
            district_kms[district] = already_walked_in_district + dist
            continue
        end

        # add distance from part[1] to next_node to the district of next node 
        next_node = nodes[start_id+1]
        first_point = get_gps_point(nodes, part[1], trans, rev_trans)
        dist = euclidean_distance(first_point, LLA(next_node.lat, next_node.lon)) / 1000
        internal_node_id = city_map.osm_id_to_node_id[next_node.id]
        district = city_map.nodes_to_district_name[internal_node_id]
        already_walked_in_district = get(district_kms, district, 0.0)
        district_kms[district] = already_walked_in_district + dist

        # add distance from part[2] to prev_idx to the district of the previous node 
        prev_node = nodes[stop_id-1]
        last_point = get_gps_point(nodes, part[2], trans, rev_trans)
        dist = euclidean_distance(last_point, LLA(prev_node.lat, prev_node.lon)) / 1000
        internal_node_id = city_map.osm_id_to_node_id[prev_node.id]
        district = city_map.nodes_to_district_name[internal_node_id]
        already_walked_in_district = get(district_kms, district, 0.0)
        district_kms[district] = already_walked_in_district + dist

        for i in start_id+1:stop_id-1
            node = nodes[i]
            next_node = nodes[i+1]
            internal_node_id = city_map.osm_id_to_node_id[node.id]
            internal_next_node_id = city_map.osm_id_to_node_id[next_node.id]
            district = city_map.nodes_to_district_name[internal_node_id]
            next_district = city_map.nodes_to_district_name[internal_next_node_id]
            dist = euclidean_distance(LLA(node.lat, node.lon), LLA(next_node.lat, next_node.lon)) / 1000
            already_walked_in_district = get(district_kms, district, 0.0)
            district_kms[district] = already_walked_in_district + dist
            if district != next_district
                already_walked_in_next_district = get(district_kms, next_district, 0.0)
                district_kms[next_district] = already_walked_in_next_district + dist
            end
        end
    end
    return district_kms
end

function get_missing_district_parts(city_map::AbstractSimpleMap, walked_parts::WalkedParts, district_name::Symbol)
    street_names = Set{String}()
    for way in city_map.ways
        check_way = false
        for node in way.nodes 
            internal_node_id = city_map.osm_id_to_node_id[node.id]
            district = city_map.nodes_to_district_name[internal_node_id]
            if district == district_name 
                check_way = true 
                break
            end
        end
        check_way || continue
        
        walked_dist, complete_dist = get_walked_way(walked_parts, way)
        if walked_dist != complete_dist
            @show way.name 
            @show walked_dist/complete_dist
            push!(street_names, way.name)
        end
    end
    @show street_names
end

function get_district_kms(city_map::AbstractSimpleMap, walked_ways::Vector{WalkedWay}; filter_fct=(way)->EverySingleStreet.iswalkable_road(way))
    district_kms = OrderedDict{Symbol, Float64}()
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)
    for walked_way in walked_ways
        filter_fct(walked_way.way) || continue
        update_district_walked!(district_kms, city_map, walked_way, trans, rev_trans)
    end
    sort!(district_kms; byvalue=true, rev=true)
    return district_kms
end

function get_district_kms(city_map::AbstractSimpleMap; filter_fct=(way)->EverySingleStreet.iswalkable_road(way))
    district_kms = OrderedDict{Symbol, Float64}()
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)
    for way in city_map.ways
        filter_fct(way) || continue
        walked_way = WalkedWay(way, [(0.0, total_length(way))])
        update_district_walked!(district_kms, city_map, walked_way, trans, rev_trans)
    end
    sort!(district_kms; byvalue=true, rev=true)
    return district_kms
end

function get_walked_district_perc(city_map::AbstractSimpleMap, walked_ways::Vector{WalkedWay})
    walked_district_kms = get_district_kms(city_map, walked_ways)
    district_kms = get_district_kms(city_map)
    district_perc = OrderedDict{Symbol, Float64}()
    for district in keys(district_kms)
        if !haskey(walked_district_kms, district)
            district_perc[district] = 0.0
            continue 
        end
        district_perc[district] = 100 * (walked_district_kms[district] / district_kms[district])
    end

    sort!(district_perc; byvalue=true, rev=true)
    return district_perc
end
