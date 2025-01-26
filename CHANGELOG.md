## Unreleased
### Bugfixes
- In certain situations the gps points might be off 
- The too fast check is only applied if the distance is over 50m
    - For shorter ones the GPS deviation plays a bigger role and speed calculation might be off
- Free memory whenever `LightXML` is used
- Shortest path has non walkable roads is now fixed with earlier by using shortest path only walkable calculation
    - This fixes a bug for which splitting still resulted in shortest path via non walkable
### Other
- Fill some small gaps

## v0.1.2 (31st of December 2024)
### Bugfixes
- District kms are calculated more correctly
    - Still a section can be counted to two districts if one node is inside one and the next inside another district
- If a point has no candidates but the point before or after has one:
    Try to find candidates in between to cover as much walked road as possible [#64](https://github.com/Wikunia/EverySingleStreet.jl/pull/64)
### Other
- Add short connectors [#60](https://github.com/Wikunia/EverySingleStreet.jl/issues/60)

## v0.1.1 (28th of November 2024)
### Bugfixes: 
- shortest path contained non walkable roads but walkable roads were feasible 
### Other
- Extend walked ways to fill GPS error gaps
- Extend roundabouts

## v0.1.0 (12th of July 2024)
- Initial version with basic functionality like:
    - map matching
    - statistics on districts