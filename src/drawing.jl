"""
    draw(walked_parts::WalkedParts, gps_points::Vector{GPSPoint}, fname; color="black", gps_opacity=0.4, line_width=7)

Draw the given walked parts on a transparent background as well as the path walked (given the gps_points) as well.
One can define the color as well as the opacity for the gps path and the overall line width.
"""
function draw(walked_parts::WalkedParts, gps_points::Vector{GPSPoint}, fname; color="black", gps_opacity=0.4, line_width=7)
    centroid = get_centroid(walked_parts)
    trans = ENUfromLLA(centroid, wgs84)
    ways = get_ways_from_walked_parts(walked_parts)
    paths = Vector{Vector{Point2{Float32}}}()
    for way in ways
        points = Vector{Point2{Float32}}() 
        for node in way.nodes 
            x, y = getxy_from_lat_lon(node.lat, node.lon, trans)
            push!(points, Point2{Float32}(x,y)) 
        end
        push!(paths, points)
    end
    gps_points = simplify(gps_points, 5)
    gps_path = Vector{Point2{Float32}}()
    for gps_point in gps_points
        x, y = getxy_from_lat_lon(gps_point.pos.lat, gps_point.pos.lon, trans)
        push!(gps_path, Point2{Float32}(x,y)) 
    end

    min_x, max_x = extrema(p->p[1], gps_path)
    min_y, max_y = extrema(p->p[2], gps_path)
    width = max(abs(max_x), abs(min_x))*2
    height = max(abs(max_y), abs(min_y))*2
    img_dim = 700

    factor = max(width/img_dim, height/img_dim)
    Luxor.Drawing(800, 800, fname)
    Luxor.origin()
    Luxor.sethue(color)
    Luxor.scale(1/factor, -1/factor)
    Luxor.setline(line_width)
    Luxor.setopacity(gps_opacity)
    lpoints = [Luxor.Point(p...) for p in gps_path]
    Luxor.poly(lpoints, :stroke)  
    Luxor.setopacity(1)
    for path in paths 
        lpoints = [Luxor.Point(p...) for p in path]
        Luxor.poly(lpoints, :stroke)  
    end
    Luxor.finish()
    return nothing
end
