function draw_way(way::Way, trans)
    setopacity(0.5)
    curve = [Point(getxy_from_lat_lon(node.lat, node.lon, trans)...) for node in way.nodes]
    poly(curve, :stroke)
end

function draw_path(path, trans)
    curve = [Point(getxy_from_lat_lon(node.lat, node.lon, trans)...) for node in path]
    poly(curve, :stroke)
end

function draw_path(candidates::Vector{Candidate}, trans)
    draw_path([c.lla for c in candidates], trans)
end

function draw_path(streetpath::StreetPath, trans)
    for segment in streetpath.segments
        draw_segment(segment, trans)
    end
end

function draw_line(lla1, lla2, trans)
    p1 = Point(getxy_from_lat_lon(lla1.lat, lla1.lon, trans))
    p2 = Point(getxy_from_lat_lon(lla2.lat, lla2.lon, trans))
    line(p1, p2, :stroke)
end

function draw_segment(segment::StreetSegment, trans)
    prev_from_idx = prev_idx(segment.from)+1
    prev_to_idx = prev_idx(segment.to)
    if segment.from == segment.to 
        Luxor.circle(Point(getxy_from_lat_lon(segment.from.lla.lat, segment.from.lla.lon, trans)), 5, :fill)
    end
    if prev_from_idx > prev_to_idx
        return draw_line(segment.from.lla, segment.to.lla, trans)
    end
    path = [segment.from.lla]
    nodes = segment.from.way.nodes
    if segment.from.way_is_reverse
        nodes = reverse(nodes)
    end

    for i in prev_from_idx:prev_to_idx
        node = nodes[i]
        push!(path, LLA(node.lat, node.lon))
    end
    push!(path, segment.to.lla)
    draw_path(path, trans)
end

function prev_idx(candidate)
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    dists = [euclidean_distance(LLA(n1.lat, n1.lon), LLA(n2.lat, n2.lon)) for (n1,n2) in zip(nodes[1:end-1], nodes[2:end])]
    cum_dists = cumsum(dists)
    pos = findfirst(>=(candidate.λ), cum_dists)
    if pos === nothing 
        return length(cum_dists)
    end
    return pos
end

function draw_candidates(city_map, candidates, outpath, scale_factor=0.1)
    origin_lla = get_centroid(city_map.nodes)
    origin_lla = candidates[1].lla
    trans = ENUfromLLA(origin_lla, wgs84)
    Drawing(1920, 1080, outpath)
    origin()
    background("white")
    Luxor.scale(scale_factor,-scale_factor)
    sethue("black")
    for way in city_map.ways
        sethue("black")
        draw_way(way, trans)
    end
    sethue("green")
    for c in candidates
        p = Point(getxy_from_lat_lon(c.lla.lat, c.lla.lon, trans))
        circle(p, 10, :fill)
        # draw_way(c.way, trans)
    end
    finish()
end

function draw_map(city_map, paths, outpath; streetpaths=Vector{StreetPath}(), scale_factor=0.1)
    origin_lla = get_centroid(city_map.nodes)
    if !isempty(streetpaths)
        origin_lla = streetpaths[1].segments[1].from.lla
    end
    trans = ENUfromLLA(origin_lla, wgs84)
    Drawing(1920, 1080, outpath)
    origin()
    background("white")
    Luxor.scale(scale_factor,-scale_factor)
    sethue("black")
    for way in city_map.ways
        sethue("black")
        draw_way(way, trans)
    end
    #=
    sethue("blue")
    setline(3)
    for path in paths
        draw_path(path, trans)
    end
    =#
    sethue("green")
    setline(3)
    for streetpath in streetpaths
        draw_path(streetpath, trans)
    end
    

    finish()
end