## Unreleased
### Bugfixes
- District kms are calculated more correctly
    - Still a section can be counted to two districts if one node is inside one and the next inside another district
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