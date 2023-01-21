function draw_way(way::Way, trans)
    setopacity(0.5)
    curve = [Point(getxy_from_lat_lon(node.lat, node.lon, trans)...) for node in way.nodes]
    poly(curve, :stroke)
end

function draw_path(path, trans)
    curve = [Point(getxy_from_lat_lon(node.lat, node.lon, trans)...) for node in path]
    poly(curve, :stroke)
end


function draw_map(map, paths, outpath; way_ids=Set(), scale_factor=0.1)
    origin_lla = get_centroid(map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    Drawing(1920, 1080, outpath)
    origin()
    background("white")
    scale(scale_factor,-scale_factor)
    sethue("black")
    for way in map.ways
        if way.id in way_ids
            sethue("green")
        else 
            sethue("black")
        end
        draw_way(way, trans)
    end
    sethue("blue")
    setline(3)
    for path in paths
        draw_path(path, trans)
    end

    finish()
end