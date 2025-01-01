var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = EverySingleStreet","category":"page"},{"location":"#EverySingleStreet","page":"Home","title":"EverySingleStreet","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for EverySingleStreet.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [EverySingleStreet]","category":"page"},{"location":"#EverySingleStreet.AbstractSimpleMap","page":"Home","title":"EverySingleStreet.AbstractSimpleMap","text":"AbstractSimpleMap\n\nAbstract supertype for different map representations.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.BoundedAllShortestPaths","page":"Home","title":"EverySingleStreet.BoundedAllShortestPaths","text":"BoundedAllShortestPaths\n\nHolds pre-computed shortest paths between nodes in a city.\n\nFields\n\ng::StaticGraphs.StaticDiGraph{Int32, Int32}: The underlying StaticDiGraph{Int32, Int32} representing the street network.\ndist_mat::SparseArrays.SparseMatrixCSC{Float64, Int32}: A SparseMatrixCSC{Float64, Int32} containing the distance matrix for the nodes in a sparse matrix.\ndistances::Vector{Dict{Int32, Float64}}: A vector of dictionaries, where each dictionary maps a destination node ID to its shortest distance from the source node.\nparents::Vector{Dict{Int32, Int32}}: A vector of dictionaries, where each dictionary maps a destination node ID to its parent node on the shortest path from the source node.\ndistance::Float64: The distance for which all the shortest paths are computed\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.Candidate","page":"Home","title":"EverySingleStreet.Candidate","text":"Candidate\n\nRepresents a candidate point potentially lying on a walkable street segment.\n\nFields\n\nmeasured_point::GPSPoint: The GPS point corresponding to the candidate.\nlla::LLA: The LLA object containing the candidate's latitude and longitude.\nway::Way: The Way object representing the street segment the candidate is on\nway_is_reverse::Bool: Indicates if the candidate is positioned along the way in the reverse direction.\ndist::Float64: The distance between the measured GPS point and the candidate point on the way.\nλ::Float64: The distance along the way on which the lla sits basically the interpolation value for the way.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.District","page":"Home","title":"EverySingleStreet.District","text":"District\n\nRepresents a district within the city.\n\nFields\n\nname::Symbol: Symbolic name of the district.\npolygons::Vector{HolePolygon}: A vector of HolePolygon objects defining the geographical area of the district.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.GPSPoint","page":"Home","title":"EverySingleStreet.GPSPoint","text":"GPSPoint\n\nRepresents a single point with GPS coordinates and associated time.\n\nFields\n\npos::LLA: An LLA object containing latitude and longitude coordinates as well as altitude (not used though)\ntime::ZonedDateTime: The time the GPS point was recorded.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.GPXFile","page":"Home","title":"EverySingleStreet.GPXFile","text":"GPXFile\n\nRepresents a GPX file containing a collection of GPS points.\n\nFields:\n\nname::String: The name of the GPX file.\ngps_points::Vector{GPSPoint}: A vector of GPSPoint objects representing the GPS points stored in the file.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.HolePolygon","page":"Home","title":"EverySingleStreet.HolePolygon","text":"HolePolygon\n\nRepresents a polygon with holes used for defining districts.\n\nFields\n\nouter::Vector{Point2{Float32}}: A vector of Point2{Float32} objects defining the outer boundary of the polygon.\nholes::Vector{Vector{Point2{Float32}}}: A vector of vectors of Point2{Float32} objects defining any holes within the polygon.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.Map","page":"Home","title":"EverySingleStreet.Map","text":"Map <: AbstractSimpleMap\n\nRepresents a map containing detailed street network information and potentially pre-computed shortest paths.\n\nFields\n\nno_graph_map::NoGraphMap: A NoGraphMap instance containing core street network data.\ngraph::Union{Missing, OSMGraph}): An optional OSMGraph instance representing the underlying street network graph (may be missing if shortest paths are not pre-computed).\nbounded_shortest_paths::Union{Missing, BoundedAllShortestPaths}: An optional BoundedAllShortestPaths instance containing pre-computed shortest paths within a specific area (may be missing if not pre-computed).\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.NoGraphMap","page":"Home","title":"EverySingleStreet.NoGraphMap","text":"NoGraphMap <: AbstractSimpleMap\n\nRepresents a map where the network graph is not explicitly stored.\n\nFields\n\nosm_id_to_node_id::Dict{Int, Int}: A dictionary mapping OpenStreetMap node IDs to internal node IDs.\nosm_id_to_edge_id::Dict{Int, Int}: A dictionary mapping OpenStreetMap way (street) IDs to internal edge IDs.\nnodes_to_district_name::Vector{Symbol}: A vector associating each node with the district it belongs to (represented by a symbol).\nnodes::Vector{Node}: A vector of Node objects representing all nodes in the network.\nways::Vector{Way}: A vector of Way objects representing all ways (streets) in the network.\nwalkable_road_nodes::Vector{Bool}: A vector indicating whether each node is on a walkable road.\nosm_node_id_to_edge_ids::Dict{Int, Vector{Int}}: A dictionary mapping OpenStreetMap node IDs to a vector of edge IDs it connects to.\ndistricts::Vector{District}: A vector of District objects defining districts within the city.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.Node","page":"Home","title":"EverySingleStreet.Node","text":"Node\n\nRepresents a single node in the street network.\n\nFields\n\nid::Int: Unique identifier for the node.\nlat::Float64: Latitude coordinate of the node.\nlon::Float64`: Longitude coordinate of the node.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.StreetPath","page":"Home","title":"EverySingleStreet.StreetPath","text":"StreetPath\n\nRepresents a path through the street network defined by a sequence of connected street segments.\n\nFields\n\nname::String: A name for the path (e.g., user-defined name).\nsubpath_id::Int: Unique identifier for the path within a larger context (e.g., route).\nsegments::Vector{StreetSegment}: A vector of StreetSegment objects representing the connected street segments forming the path.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.StreetSegment","page":"Home","title":"EverySingleStreet.StreetSegment","text":"struct StreetSegment\n\nA street segment contains of two candidates which have the following property: They are both on the same way and have the same direction.\n\nFields\n\nfrom::Candidate the start Candidate of the segment.\nto::Candidate the start Candidate of the segment.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.WalkedParts","page":"Home","title":"EverySingleStreet.WalkedParts","text":"WalkedParts\n\nRepresents a collection of walked way segments for a specific area, potentially with additional information for each way.\n\nFields:\n\nnames::Dict{String, Vector{Int}}: A dictionary where keys are way names (from Way.name) and values are vectors of integers referencing corresponding WalkedWay objects in the ways field.\nways::Dict{Int, WalkedWay}: A dictionary where keys are unique identifiers and values are WalkedWay objects representing the walked way segments.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.WalkedWay","page":"Home","title":"EverySingleStreet.WalkedWay","text":"WalkedWay\n\nRepresents a way (street) segment that has been walked along, potentially with additional information about the walked parts.\n\nFields\n\nway::Way: The Way object representing the underlying street segment.\nparts::Vector{Tuple{Float64, Float64}}: A vector of tuples where each tuple represents a portion of the way that was walked (start and end distance along the way in meters).\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.Way","page":"Home","title":"EverySingleStreet.Way","text":"Way\n\nRepresents a single way (street) in the street network. Most data is coming from openstreetmap.org so for more information about the meaning of highway as an example their documentation should be checked.\n\nFields\n\nid::Int: Unique identifier for the way.\nnodes::Vector{Node}: A vector of Node objects defining the way.\nname::String: Name of the way (e.g., street name).\nhighway (String): Highway type of the way (e.g., \"motorway\", \"residential\").\nfoot (String): Access for pedestrians (e.g., \"yes\", \"no\").\naccess (String): General access information for the way.\nmeters::Float64: Total length of the way in meters.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.bounded_all_shortest_paths-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.bounded_all_shortest_paths","text":"bounded_all_shortest_paths(osm_graph, distance, nodeid_to_local, is_walkable_road)\n\nCall bounded_dijkstra for all nodes and return a BoundedAllShortestPaths object.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.bounded_dijkstra-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.bounded_dijkstra","text":"bounded_dijkstra(graph, dist_mat, source, distance)\n\nCompute the shortest distance and the parents on how to get there from the specified source up to distance away. Return distances as a dictionary pointing from destination to shortest distance as well as parents point from destination on how to get there.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.calculate_streetpath-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.calculate_streetpath","text":"calculate_streetpath(name, subpath_id, candidates, city_map)\n\nGenerate a StreetPath from a list of candidates obtained by map_path.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.calculate_walked_parts","page":"Home","title":"EverySingleStreet.calculate_walked_parts","text":"calculate_walked_parts(streetpaths::Vector{StreetPath}, city_map::AbstractSimpleMap, walked_ways=Dict{Int, WalkedWay}())\n\nReturn WalkedParts given the streetpath segments and the possible ways of the city. Can be added to already existing walked_ways. Filters everything which isn't a walkable_road inside city_ways.    Calls several function to get a better approximation of what was actually walked like      - add extra buffer at start and end of streets to not miss the last few meters of a dead end street as an example extend_walked_parts!     - closes circles like roundabouts, parts at ends of some dead end streets\n\n\n\n\n\n","category":"function"},{"location":"#EverySingleStreet.check_if_visited_node-Tuple{Any, EverySingleStreet.Node, Any}","page":"Home","title":"EverySingleStreet.check_if_visited_node","text":"check_if_visited_node(city_map, node::Node, walked_parts)\n\nReturn whether the node was visited already.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.douglas_peucker-Tuple{Vector, Any, Any, Any}","page":"Home","title":"EverySingleStreet.douglas_peucker","text":"Use a non-recursive Douglas-Peucker algorithm to simplify a polygon. Used by simplify().\n\ndouglas_peucker(pointlist::Array, start_index, last_index, epsilon)\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.download-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.download","text":"download(place_name, filepath)\n\nDownload the road network of the given place and write it as a json file to the given filepath\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.draw-Tuple{EverySingleStreet.WalkedParts, Vector{EverySingleStreet.GPSPoint}, Any}","page":"Home","title":"EverySingleStreet.draw","text":"draw(walked_parts::WalkedParts, gps_points::Vector{GPSPoint}, fname; color=\"black\", gps_opacity=0.4, line_width=7)\n\nDraw the given walked parts on a transparent background as well as the path walked (given the gps_points) as well. One can define the color as well as the opacity for the gps path and the overall line width.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.extend_walked_parts_connectors!-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.extend_walked_parts_connectors!","text":"extend_walked_parts_connectors!(walked_parts, city_map)\n\nFind all the city ways that are shorter than 2*EXTENDWALKEDWAYUPTO\nCheck if the end of those ways is part of some already walked \n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.extend_walked_parts_cycle!-Tuple{Any}","page":"Home","title":"EverySingleStreet.extend_walked_parts_cycle!","text":"extend_walked_parts_cycle!(walked_parts)\n\nExtend the walked parts by finishing a cycle if more than MIN_FILL_CYCLE_PERC perc of a cycle of length maximum MAX_FILL_CYCLE_LENGTH is walked.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.extend_walked_parts_simple!-Tuple{Any}","page":"Home","title":"EverySingleStreet.extend_walked_parts_simple!","text":"extend_walked_parts_simple!(walked_parts)\n\nExtend the walked parts by adding up to EXTEND_WALKED_WAY_UP_TO to the start and end of a walked way. As well as to fill gaps in between parts.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.filter_candidates!-Tuple{Any}","page":"Home","title":"EverySingleStreet.filter_candidates!","text":"filter_candidates!(candidates)\n\nFilter out candidates which are further away than closer_dist if there is at least candidate which is closer than closer_dist.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.find_in_between_candidates-Tuple{EverySingleStreet.AbstractSimpleMap, NearestNeighbors.KDTree, Any, EverySingleStreet.GPSPoint, EverySingleStreet.GPSPoint, Geodesy.LLA, Any}","page":"Home","title":"EverySingleStreet.find_in_between_candidates","text":"find_in_between_candidates(city_map::AbstractSimpleMap, way_tree::KDTree, id_to_way_id, p1::GPSPoint, p2::GPSPoint, origin_lla::LLA, radius)\n\nReturn candidates that are as close as possible to p2 on the line from p1 to p2.  This should be called when there are no candidates for p2 but there are for p1.  In that case we want to find candidates for the last point that we still have candidates for.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidate_on_way-Tuple{Any, Any, EverySingleStreet.Way, Any, Any}","page":"Home","title":"EverySingleStreet.get_candidate_on_way","text":"get_candidate_on_way(city_map, p, way::Way, trans, rev_trans; rev=false)\n\nGet the best candidate of point p on the given way. trans and rev_trans are transformations mapping from LLA to x,y and back. rev can be set to true to reverse the direction of way.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidate_probability-Tuple{EverySingleStreet.Candidate}","page":"Home","title":"EverySingleStreet.get_candidate_probability","text":"get_candidate_probability(candidate::Candidate)\n\nReturn the emission probability of the given candidate and the standard deviation of the gps error.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidates-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.get_candidates","text":"get_candidates(map, path)\n\nGet a list of Candidate for each gps point in path given the underlying city_map\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidates-Tuple{EverySingleStreet.AbstractSimpleMap, NearestNeighbors.KDTree, Any, EverySingleStreet.GPSPoint, Geodesy.LLA, Any}","page":"Home","title":"EverySingleStreet.get_candidates","text":"get_candidates(city_map::AbstractSimpleMap, way_tree::KDTree, id_to_way_id, p::GPSPoint, origin_lla::LLA, radius)\n\nGet canddiates for a single point given already a kdtree which can map the point to possible candidate ways and so on. This function is used by find_in_between_candidates. \n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidates_from_idx-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.get_candidates_from_idx","text":"get_candidates_from_idx(vec_candidates, candidate_idxs)\n\nReturn candidates from an nested vector of candiates like [[c1,c2],[c3]] and candiate_idxs which are 1d like [1,3] would return [c1, c3].\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_districts-Tuple{Nothing}","page":"Home","title":"EverySingleStreet.get_districts","text":"get_districts(geojson_fpath)\n\nParse the districts in the given file which needs to have :geometry and :Stadtteil as properties. Return a vector of District\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_first_way_segment-Tuple{Any, EverySingleStreet.Map}","page":"Home","title":"EverySingleStreet.get_first_way_segment","text":"get_first_way_segment(sp, city_map::Map)\n\nReturn a way that includes the longest first part of the given shortest path  as well as the remaining shortest path which isn't part of this way. The return format consists of:\n\nA named tuple describing the best way for the first segment: (way::Way, rev::Bool, from::Int, to::Int)\nThe remaining shortest path nodes\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_interpolation_point-Tuple{Any, Any, Float64}","page":"Home","title":"EverySingleStreet.get_interpolation_point","text":"get_interpolation_point(a, b, t::Float64)\n\nGet the value of a+t*(b-a)\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_interpolation_val-Tuple{GeometryBasics.Point2, GeometryBasics.Point2, GeometryBasics.Point2}","page":"Home","title":"EverySingleStreet.get_interpolation_val","text":"get_interpolation_val(p::Point2, a::Point2, b::Point2)\n\nGet an interpolation value of point p between a line segment from a to b where a+t*(b-a) describes the point closes to p. Return p which is between 0 and 1.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_lla-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.get_lla","text":"get_lla(p, trans)\n\nGet a LLA formatted position from a point and a Geodesy transformation.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_matching_candidates-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.get_matching_candidates","text":"get_matching_candidates(city_map, way_ids, p, origin_lla)\n\nGet matching candidates for a given point p and a list of possible way ids. Only return candidates which have a maximum_dist (in m) to the way.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_next_node_id-Tuple{EverySingleStreet.Candidate}","page":"Home","title":"EverySingleStreet.get_next_node_id","text":"get_next_node_id(candidate::Candidate)\n\nGet the next node id given a candidate. This depends on way_is_reverse of the candidate. Can be used to create a StreetSegment.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_prev_node_id-Tuple{EverySingleStreet.Candidate}","page":"Home","title":"EverySingleStreet.get_prev_node_id","text":"get_prev_node_id(candidate::Candidate)\n\nGet the previous node id given a candidate. This depends on way_is_reverse of the candidate. Can be used to create a StreetSegment.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_segments-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.get_segments","text":"get_segments(city_map, current_candidate, next_candidate, sp)\n\nGiven two candidates which can't be directly connected and a list of ids which form the shortest path between those candidates this function returns a list of StreetSegment which connect the two candidates.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_shortest_path-Tuple{EverySingleStreet.BoundedAllShortestPaths, Any, Any}","page":"Home","title":"EverySingleStreet.get_shortest_path","text":"get_shortest_path(bounded_shortest_paths::BoundedAllShortestPaths, from, to)\n\nReturn the same output as a_star(g, from, to, dist_mat) but use the cache from BoundedAllShortestPaths\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_shortest_path-Tuple{EverySingleStreet.Map, Any, Any}","page":"Home","title":"EverySingleStreet.get_shortest_path","text":"get_shortest_path(city_map::Map, sp_from_id, sp_to_id)\n\nReturn the shortest path from from to to. Has the same output as shortest_path from the LightOSM package  but uses the BoundedAllShortestPaths cache when it exists.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_way_kdtree-Tuple{Any}","page":"Home","title":"EverySingleStreet.get_way_kdtree","text":"get_way_kdtree(city_map)\n\nCreate a KDTree for all ways in the given city map. Return\n\na mapping from id to the way id\nthe kd tree\nThe radius for a reasonable in range search (51% of the longest line between two nodes)\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.getxy_from_lat_lon-Tuple{Any, Any, Any}","page":"Home","title":"EverySingleStreet.getxy_from_lat_lon","text":"getxy_from_lat_lon(lat, lon, trans)\n\nReturn x, y position coordinates from lat and lon given a Geodesy trensformation.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.isfullywalked-Tuple{Any, EverySingleStreet.Way}","page":"Home","title":"EverySingleStreet.isfullywalked","text":"isfullywalked(walked_parts, way::Way)\n\nReturn whether the way is already completely walked.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.map_nodes_to_district-Tuple{Vector{EverySingleStreet.Node}, Any}","page":"Home","title":"EverySingleStreet.map_nodes_to_district","text":"map_nodes_to_district(nodes::Vector{Node}, geojson_fpath)\n\nReturn a vector of district names (as Symbol) that points each Node to its respective district\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.map_path-Tuple{Any, EverySingleStreet.GPXFile}","page":"Home","title":"EverySingleStreet.map_path","text":"map_path(city_map, gpxfile)\n\nMap a path of gps points to the best matching candidates for the path.\n\nFilter out points on the path which are closer together than 25m\nCompute candidates for each point in the remaining path\nCompute emission probabilties\nCompute transition_probabilities\nUse the Viterbi algorithm to compute the most likely path of the candidates\n\nReturn a list of lists of candidates as a path might not continously have candidates. Then it's split up into several parts.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.midpoint-Tuple{EverySingleStreet.GPSPoint, EverySingleStreet.GPSPoint}","page":"Home","title":"EverySingleStreet.midpoint","text":"midpoint(p1::GPSPoint, p2::GPSPoint)\n\nCalculates the midpoint between two GPS points.\n\nThis function takes two GPSPoint objects as input, each containing a position represented by LLA (Latitude, Longitude, Altitude) coordinates and a timestamp represented by ZonedDateTime. It calculates the midpoint in ECEF (Earth-Centered, Earth-Fixed) coordinates for accuracy, then converts the result back to LLA. The time of the midpoint is also interpolated linearly.\n\nArguments\n\np1::GPSPoint: The first GPS point.\np2::GPSPoint: The second GPS point.\n\nReturns\n\nA GPSPoint object representing the midpoint between p1 and p2.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.node_on_way_to_candidate-Tuple{Int64, EverySingleStreet.Way}","page":"Home","title":"EverySingleStreet.node_on_way_to_candidate","text":"node_on_way_to_candidate(idx::Int, way::Way; rev=false)\nnode_on_way_to_candidate(node::Node, way::Way; rev=false)\n\nGiven either the idx of a node on a way like way.nodes[idx] if not reversed or the node itself return the matching candidate by computing λ the time of the candidate is set to now in UTC.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.parse_map","page":"Home","title":"EverySingleStreet.parse_map","text":"parse_map(fpath)\n\nReturn a Map object from the given json file path that was created using the download function.\n\n\n\n\n\n","category":"function"},{"location":"#EverySingleStreet.point_linesegment_distance-Tuple{GeometryBasics.Point2, GeometryBasics.Point2, GeometryBasics.Point2}","page":"Home","title":"EverySingleStreet.point_linesegment_distance","text":"point_linesegment_distance(p::Point2, a::Point2, b::Point2)\n\nGet the smallest distance between point p and the line segment from a to b.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.pointlinedistance-Tuple{GeometryBasics.Point2, GeometryBasics.Point2, GeometryBasics.Point2}","page":"Home","title":"EverySingleStreet.pointlinedistance","text":"pointlinedistance(p::Point2, a::Point2, b::Point2)\n\nFind the distance between a point p and a line between two points a and b.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.shortest_candidate_path-Tuple{EverySingleStreet.Candidate, EverySingleStreet.Candidate, Any}","page":"Home","title":"EverySingleStreet.shortest_candidate_path","text":"shortest_candidate_path(from::Candidate, to::Candidate, city_map)\n\nReturn the shortest path as node ids from one candidate to another given a city map.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.shortest_path_distance-Tuple{EverySingleStreet.Candidate, EverySingleStreet.Candidate, Any}","page":"Home","title":"EverySingleStreet.shortest_path_distance","text":"shortest_path_distance(from::Candidate, to::Candidate, city_map)\n\nReturn the shortest path distance based on equation (9) in the fmm paper\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.simplify","page":"Home","title":"EverySingleStreet.simplify","text":"Simplify a polygon:\n\nsimplify(pointlist::Array, detail=0.1)\n\ndetail is the maximum approximation error of simplified polygon.\n\n\n\n\n\n","category":"function"},{"location":"#EverySingleStreet.transition_probability-Tuple{EverySingleStreet.Candidate, EverySingleStreet.Candidate, Any}","page":"Home","title":"EverySingleStreet.transition_probability","text":"transition_probability(from::Candidate, to::Candidate, city_map)\n\nReturn the transitionprobability based on equation 10 in the  [fmm paper](https://www.tandfonline.com/doi/pdf/10.1080/13658816.2017.1400548?casatoken=RrxNeRXfRfkAAAAA:2IA6z4Pu-tKEcSaK44AUgQxDc-XUCBs8CSZI1qGNbUj6CpUMyA8suDUpnZ1WO3lHEUFuk1lk3s4wJtM)\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.update_district_walked!-Tuple{AbstractDict{Symbol, Float64}, EverySingleStreet.AbstractSimpleMap, EverySingleStreet.WalkedWay, Any, Any}","page":"Home","title":"EverySingleStreet.update_district_walked!","text":"update_district_walked!(district_kms::AbstractDict{Symbol, Float64}, city_map::AbstractSimpleMap, walkedway::WalkedWay)\n\nUpdate the total walked distance of each district after a WalkedWay is added to the map.\n\nArguments\n\ndistrict_kms::AbstractDict{Symbol, Float64}: A dictionary mapping district names (as Symbol) to their total walked distances in kilometers.\ncity_map::AbstractSimpleMap: The city map object that stores nodes and edges.\nwalkedway::WalkedWay: The walked way object added to the map.\n\n\n\n\n\n","category":"method"}]
}