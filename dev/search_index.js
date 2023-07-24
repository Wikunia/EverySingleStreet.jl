var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = EverySingleStreet","category":"page"},{"location":"#EverySingleStreet","page":"Home","title":"EverySingleStreet","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for EverySingleStreet.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [EverySingleStreet]","category":"page"},{"location":"#EverySingleStreet.StreetSegment","page":"Home","title":"EverySingleStreet.StreetSegment","text":"struct StreetSegment\n\nA street segment contains of two candidates which have the following property: They are both on the same way and have the same direction.\n\n\n\n\n\n","category":"type"},{"location":"#EverySingleStreet.bounded_all_shortest_paths-Tuple{Any, Any, Any}","page":"Home","title":"EverySingleStreet.bounded_all_shortest_paths","text":"bounded_all_shortest_paths(g, dist_mat, distance)\n\nCall bounded_dijkstra for all nodes and return a BoundedAllShortestPaths object.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.bounded_dijkstra-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.bounded_dijkstra","text":"bounded_dijkstra(graph, dist_mat, source, distance)\n\nCompute the shortest distance and the parents on how to get there from the specified source up to distance away. Return distances as a dictionary pointing from destination to shortest distance as well as parents point from destination on how to get there.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.calculate_streetpath-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.calculate_streetpath","text":"calculate_streetpath(ame, subpath_id, candidates, city_map)\n\nGenerate a StreetPath from a list of candidates obtained by map_path.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.download-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.download","text":"download(place_name, filepath)\n\nDownload the road network of the given place and write it as a json file to the given filepath\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.filter_candidates!-Tuple{Any}","page":"Home","title":"EverySingleStreet.filter_candidates!","text":"filter_candidates!(candidates; closer_dist=25)\n\nFilter out candidates which are further away than closer_dist if there is at least candidate which is closer than closer_dist.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.filter_path-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.filter_path","text":"filter_path(gps_points, dist)\n\nReturn a new path for which the minimum distance between two consecutive points is at least dist.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidate_on_way-Tuple{Any, Any, EverySingleStreet.Way, Any, Any}","page":"Home","title":"EverySingleStreet.get_candidate_on_way","text":"get_candidate_on_way(city_map, p, way::Way, trans, rev_trans; rev=false)\n\nGet the best candidate of point p on the given way. trans and rev_trans are transformations mapping from LLA to x,y and back. rev can be set to true to reverse the direction of way.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidate_probability-Tuple{EverySingleStreet.Candidate}","page":"Home","title":"EverySingleStreet.get_candidate_probability","text":"get_candidate_probability(candidate::Candidate; sigma=10)\n\nReturn the emission probability of the given candidate and the standard deviation of the gps error.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidates-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.get_candidates","text":"get_candidates(map, path)\n\nGet a list of Candidate for each gps point in path given the underlying city_map\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_candidates_from_idx-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.get_candidates_from_idx","text":"get_candidates_from_idx(vec_candidates, candidate_idxs)\n\nReturn candidates from an nested vector of candiates like [[c1,c2],[c3]] and candiate_idxs which are 1d like [1,3] would return [c1, c3].\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_districts-Tuple{Any}","page":"Home","title":"EverySingleStreet.get_districts","text":"get_districts(geojson_fpath)\n\nParse the districts in the given file which needs to have :geometry and :Stadtteil as properties. Return a vector of District\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_interpolation_point-Tuple{Any, Any, Float64}","page":"Home","title":"EverySingleStreet.get_interpolation_point","text":"get_interpolation_point(a, b, t::Float64)\n\nGet the value of a+t*(b-a)\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_interpolation_val-Tuple{Luxor.Point, Luxor.Point, Luxor.Point}","page":"Home","title":"EverySingleStreet.get_interpolation_val","text":"get_interpolation_val(p::Point, a::Point, b::Point)\n\nGet an interpolation value of point p between a line segment from a to b where a+t*(b-a) describes the point closes to p. Return p which is between 0 and 1.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_lla-Tuple{Any, Any}","page":"Home","title":"EverySingleStreet.get_lla","text":"get_lla(p, trans)\n\nGet a LLA formatted position from a point and a Geodesy transformation.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_matching_candidates-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.get_matching_candidates","text":"get_matching_candidates(city_map, way_ids, p, origin_lla; maximum_dist=30)\n\nGet matching candidates for a given point p and a list of possible way ids. Only return candidates which have a maximum_dist (in m) to the way.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_next_node_id-Tuple{EverySingleStreet.Candidate}","page":"Home","title":"EverySingleStreet.get_next_node_id","text":"get_next_node_id(candidate::Candidate)\n\nGet the next node id given a candidate. This depends on way_is_reverse of the candidate. Can be used to create a StreetSegment.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_prev_node_id-Tuple{EverySingleStreet.Candidate}","page":"Home","title":"EverySingleStreet.get_prev_node_id","text":"get_prev_node_id(candidate::Candidate)\n\nGet the previous node id given a candidate. This depends on way_is_reverse of the candidate. Can be used to create a StreetSegment.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_segments-NTuple{4, Any}","page":"Home","title":"EverySingleStreet.get_segments","text":"get_segments(city_map, current_candidate, next_candidate, sp)\n\nGiven two candidates which can't be directly connected and a list of ids which form the shortest path between those candidates this function returns a list of StreetSegment which connect the two candidates.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_shortest_path-Tuple{EverySingleStreet.BoundedAllShortestPaths, Any, Any}","page":"Home","title":"EverySingleStreet.get_shortest_path","text":"get_shortest_path(bounded_shortest_paths::BoundedAllShortestPaths, from, to)\n\nReturn the same output as a_star(g, from, to, dist_mat) but use the cache from BoundedAllShortestPaths\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_shortest_path-Tuple{EverySingleStreet.Map, Any, Any}","page":"Home","title":"EverySingleStreet.get_shortest_path","text":"get_shortest_path(city_map::Map, sp_from_id, sp_to_id)\n\nReturn the shortest path from from to to. Has the same output as shortest_path from the LightOSM package  but uses the BoundedAllShortestPaths cache when it exists.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.get_way_kdtree-Tuple{Any}","page":"Home","title":"EverySingleStreet.get_way_kdtree","text":"get_way_kdtree(city_map)\n\nCreate a KDTree for all ways in the given city map. Return\n\na mapping from id to the way id\nthe kd tree\nThe radius for a reasonable in range search (51% of the longest line between two nodes)\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.getxy_from_lat_lon-Tuple{Any, Any, Any}","page":"Home","title":"EverySingleStreet.getxy_from_lat_lon","text":"getxy_from_lat_lon(lat, lon, trans)\n\nReturn x, y position coordinates from lat and lon given a Geodesy trensformation.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.map_nodes_to_district-Tuple{Vector{EverySingleStreet.Node}, Any}","page":"Home","title":"EverySingleStreet.map_nodes_to_district","text":"map_nodes_to_district(nodes::Vector{Node}, geojson_fpath)\n\nReturn a vector of district names (as Symbol) that points each Node to its respective district\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.map_path-Tuple{Any, EverySingleStreet.GPXFile}","page":"Home","title":"EverySingleStreet.map_path","text":"map_path(city_map, gpxfile)\n\nMap a path of gps points to the best matching candidates for the path.\n\nFilter out points on the path which are closer together than 25m\nCompute candidates for each point in the remaining path\nCompute emission probabilties\nCompute transition_probabilities\nUse the Viterbi algorithm to compute the most likely path of the candidates\n\nReturn a list of lists of candidates as a path might not continously have candidates. Then it's split up into several parts.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.parse_map","page":"Home","title":"EverySingleStreet.parse_map","text":"parse_map(fpath)\n\nReturn a Map object from the given json file path that was created using the download function.\n\n\n\n\n\n","category":"function"},{"location":"#EverySingleStreet.point_linesegment_distance-Tuple{Luxor.Point, Luxor.Point, Luxor.Point}","page":"Home","title":"EverySingleStreet.point_linesegment_distance","text":"point_linesegment_distance(p::Point, a::Point, b::Point)\n\nGet the smallest distance between point p and the line segment from a to b.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.shortest_candidate_path-Tuple{EverySingleStreet.Candidate, EverySingleStreet.Candidate, Any}","page":"Home","title":"EverySingleStreet.shortest_candidate_path","text":"shortest_candidate_path(from::Candidate, to::Candidate, city_map)\n\nReturn the shortest path as node ids from one candidate to another given a city map.\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.shortest_path_distance-Tuple{EverySingleStreet.Candidate, EverySingleStreet.Candidate, Any}","page":"Home","title":"EverySingleStreet.shortest_path_distance","text":"shortest_path_distance(from::Candidate, to::Candidate, city_map)\n\nReturn the shortest path distance based on equation (9) in the fmm paper\n\n\n\n\n\n","category":"method"},{"location":"#EverySingleStreet.transition_probability-Tuple{EverySingleStreet.Candidate, EverySingleStreet.Candidate, Any}","page":"Home","title":"EverySingleStreet.transition_probability","text":"transition_probability(from::Candidate, to::Candidate, city_map)\n\nReturn the transitionprobability based on equation 10 in the  [fmm paper](https://www.tandfonline.com/doi/pdf/10.1080/13658816.2017.1400548?casatoken=RrxNeRXfRfkAAAAA:2IA6z4Pu-tKEcSaK44AUgQxDc-XUCBs8CSZI1qGNbUj6CpUMyA8suDUpnZ1WO3lHEUFuk1lk3s4wJtM)\n\n\n\n\n\n","category":"method"}]
}
