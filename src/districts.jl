"""
    get_districts(geojson_fpath)

Parse the districts in the given file which needs to have :geometry and :Stadtteil as properties.
Return a vector of [`District`](@ref)
"""
function get_districts(geojson_fpath)
    districts = GeoJSON.read(read(geojson_fpath))
    df_districts = DataFrame(districts) 
    districts = Vector{District}()
    for row in eachrow(df_districts)
        hpolygons = Vector{HolePolygon}()
        geometry = row[:geometry]
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
        push!(districts, District(Symbol(row[:Stadtteil]), hpolygons))
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

function update_district_walked!(district_kms::AbstractDict{Symbol, Float64}, city_map::Map, walkedway::WalkedWay)
    way = walkedway.way
    nodes = way.nodes
    for part in walkedway.parts
        start_id = prev_idx(nodes, part[1])
        stop_id = prev_idx(nodes, part[2])+1
        for i in start_id:stop_id-1
            node = nodes[i]
            next_node = nodes[i+1]
            internal_node_id = city_map.osm_id_to_node_id[node.id]
            next_internal_node_id = city_map.osm_id_to_node_id[next_node.id]
            district = city_map.nodes_to_district_name[internal_node_id]
            next_district = city_map.nodes_to_district_name[next_internal_node_id]
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

function get_missing_district_parts(city_map::Map, walked_parts::WalkedParts, district_name::Symbol)
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

function get_district_kms(city_map::Map, walked_ways::Vector{WalkedWay})
    district_kms = OrderedDict{Symbol, Float64}()
    for walked_way in walked_ways
        update_district_walked!(district_kms, city_map, walked_way)
    end
    sort!(district_kms; byvalue=true, rev=true)
    return district_kms
end

function get_district_kms(city_map::Map)
    district_kms = OrderedDict{Symbol, Float64}()
    for way in city_map.ways
        walked_way = WalkedWay(way, [(0.0, total_length(way))])
        update_district_walked!(district_kms, city_map, walked_way)
    end
    sort!(district_kms; byvalue=true, rev=true)
    return district_kms
end

function get_walked_district_perc(city_map::Map, walked_ways::Vector{WalkedWay})
    walked_district_kms = get_district_kms(city_map, walked_ways)
    district_kms = get_district_kms(city_map)
    district_perc = OrderedDict{Symbol, Float64}()
    for district in keys(walked_district_kms)
        district_perc[district] = 100 * (walked_district_kms[district] / district_kms[district])
    end

    sort!(district_perc; byvalue=true, rev=true)
    return district_perc
end
