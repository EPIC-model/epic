## EPIC version 0.12.0
* [replace HDF5 with NetCDF](https://github.com/matt-frey/epic/pull/325)
## EPIC version 0.11.0
#### New features and enhancements
Development of the 3D model (see [#230](https://github.com/matt-frey/epic/issues/230)):
* add eigenvalue solver for real symmetric matrices (see [#229](https://github.com/matt-frey/epic/pull/229))
* add ellipsoid module (see [#233](https://github.com/matt-frey/epic/pull/233))
* add 3D parcel splitting (see [#237](https://github.com/matt-frey/epic/pull/237))
* add 3D parcel initialisation (see [#238](https://github.com/matt-frey/epic/pull/238))
* add 2D FFT (see [#240](https://github.com/matt-frey/epic/pull/240))
* add 3D parcel correction (see [#241](https://github.com/matt-frey/epic/pull/241))
* add 3D parcel merging (see [#243](https://github.com/matt-frey/epic/pull/243))
* add 3D RK4 (see [#244](https://github.com/matt-frey/epic/pull/244))
* add 3D vor2vel (see [#247](https://github.com/matt-frey/epic/pull/247))
* add 3D vorticity tendency calculation (see [#248](https://github.com/matt-frey/epic/pull/248) and [#249](https://github.com/matt-frey/epic/pull/249))
* add 3D adaptive time step estimate (see [#250](https://github.com/matt-frey/epic/pull/250))
* add changelog (see [#255](https://github.com/matt-frey/epic/pull/255))
* add unit tests (see [#242](https://github.com/matt-frey/epic/pull/242), [#246](https://github.com/matt-frey/epic/pull/246))
* [change memory layout](https://github.com/matt-frey/epic/pull/264)
* improve performance (see [#266](https://github.com/matt-frey/epic/pull/266))
* add restarting feature (see [#278](https://github.com/matt-frey/epic/pull/278))
#### Bug fixes
* [fix bug in Straka diagnostics](https://github.com/matt-frey/epic/pull/265)
* [fix calculation of get_delx and get_dely for non-centred domains](https://github.com/matt-frey/epic/pull/279)
## EPIC version 0.10.8
* [fix access of datasets + add colorcet license](https://github.com/matt-frey/epic/pull/321)
## EPIC version 0.10.7
* [fields must be reset](https://github.com/matt-frey/epic/pull/315)
* [fix peref calculation](https://github.com/matt-frey/epic/pull/317)
* [improve Python scripts](https://github.com/matt-frey/epic/pull/318)
* [change setup of parcel correction tests](https://github.com/matt-frey/epic/pull/316)
* [minor plotting changes](https://github.com/matt-frey/epic/pull/319)
## EPIC version 0.10.6
* [more buoyancy and vorticity diagnostics](https://github.com/matt-frey/epic/pull/310)
* [update python scripts](https://github.com/matt-frey/epic/pull/308)
* [do not use data to set dmin and dmax](https://github.com/matt-frey/epic/pull/311)
* [fix genspec writing spectrum](https://github.com/matt-frey/epic/pull/309)
* [fix name clash](https://github.com/matt-frey/epic/pull/312)
* [minval and maxval](https://github.com/matt-frey/epic/pull/313)
* [use heap memory](https://github.com/matt-frey/epic/pull/314)
## EPIC version 0.10.5
* [change order of multi-dimensional arrays](https://github.com/matt-frey/epic/pull/302)
* [remove merger dump](https://github.com/matt-frey/epic/pull/304)
## EPIC version 0.10.4
* [fix name of executable](https://github.com/matt-frey/epic/pull/299)
* [fix name of executable](https://github.com/matt-frey/epic/pull/301)
## EPIC version 0.10.3
* [restarting EPIC](https://github.com/matt-frey/epic/pull/278)
* [fix calculation of get_delx and get_dely for non-centred domains](https://github.com/matt-frey/epic/pull/279)
* [update pillow dependency](https://github.com/matt-frey/epic/pull/280)
* [update pyyaml dependency](https://github.com/matt-frey/epic/pull/281)
* [update urllib3 dependency](https://github.com/matt-frey/epic/pull/282)
* [fix get_delx for non-centred domains](https://github.com/matt-frey/epic/pull/283)
* [update python dependencies](https://github.com/matt-frey/epic/pull/285)
* [use proper interpolation name](https://github.com/matt-frey/epic/pull/287)
* [new directory structure in main like in branch develop](https://github.com/matt-frey/epic/pull/291)
* [combine to array](https://github.com/matt-frey/epic/pull/293)
* [add initial time like in the 3D model](https://github.com/matt-frey/epic/pull/295)
* [apply periodic BCs to parcel mergers](https://github.com/matt-frey/epic/pull/298)
## EPIC version 0.10.2
* [fix import in README](https://github.com/matt-frey/epic/commit/fa36d763f9ee9d16c60b9bca58a5bc60e5464105)
## EPIC version 0.10.1
* [use plain symbols for energy quantities in Python plots](https://github.com/matt-frey/epic/commit/9337c95b6851f7f0b7546f7ff40d3acdbad78844)
## EPIC version 0.10.0
* development of the 2D model (see [#1](https://github.com/matt-frey/epic/pull/1) to [#225](https://github.com/matt-frey/epic/pull/225))
