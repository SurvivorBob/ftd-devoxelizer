FTD Devoxelizer
===============

Intelligently converts voxel files into FTD blueprints.

![Example image](img/condor%20overview.jpg)

Prerequisites
-------------
* Python 3.8+ (not tested on older Python 3.x versions)
* A JSON voxel file to convert (described below)
* A "donor" blueprint to copy an author tag from

What
---
This command-line tool takes in a file with a solid 
[voxel representation](https://drububu.com/miscellaneous/voxelizer/?out=jso)
of some object, and produces a From the Depths blueprint with the
following properties:

* Voxels are sorted into concentric layers starting from outside, and
  the resulting blocks are colored by layer depth.
* The outermost layer is partially smoothed.
* Consecutive runs of cubes within a layer are merged into beams in
  greedy, mostly arbitrary fashion.

Below is an example of the layering behaviour, demonstrated by several
explosion craters made on the vehicle. As the layers are independently
colored, they may be easily mass-replaced in FTD through the armour
refit tool.

![Example image](img/condor%20cratered.jpg)

We copy the author tag from an existing "donor blueprint" instead of
trying to originate a tag on our own. Note that this also copies the
object ID, so you may wish to create one donor blueprint per unique
ship design you plan to create with this tool. Creating a donor
blueprint is as simple as creating a new vehicle with a single block
and saving it with a new name.

Usage
---
```
python3 convert.py [-h] [--thickness N] [--no-smoothing] [--no-beams]
                  input_file donor_blueprint output_blueprint
```

You must specify the following:
* `input_file`: a JSON file representing a __solid__ voxel 
  representation of a desired solid, containing a single JSON object
  with the following schema:
  * `dimension` is a list containing a single JSON object with the
    keys `width`, `height`, and `depth` mapping to the maximum x, y, and
    z values of voxels
  * `voxels` is a list containing one JSON object per voxel, each object
    with the keys `x`, `y` and `z` mapping to the x, y, and z positions
    of that voxel
  * +x corresponds to right, +y corresponds to forward, +z corresponds
    to up (TODO: make this configurable at runtime)
* `donor_blueprint`: a FTD blueprint file containing an authoring tag
  (author and object ID information) you want to associate with the
  final blueprint
* `output_blueprint`: a filename to store the resulting blueprint

[This page](https://drububu.com/miscellaneous/voxelizer/?out=jso) will
convert many 3D filetypes into the correct JSON format, if it still
exists on the Internet. To produce a solid voxel representation, set
the shell thickness to 2048.

You may also specify the following:
* `--thickness N` (default 2): changes the number of layers produced
  in the final blueprint. Useful for automatically defining the 
  thickness of the resulting armor. Max 32 because we run out of colors
  otherwise.
* `--no-smoothing`: if set, skips the step where we replace some blocks
  in the outermost layer with slopes and corners (if you want to do this
  differently entirely yourself)
* `--no-beams`: if set, skips the step where we opportunistically
  convert runs of single blocks into beams (if you want all single
  blocks for some reason)

Example usage
-------------
```
python convert.py "/path/to/some/solid.json" "/path/to/Constructables/some_donor.blueprint" "/path/to/Constructables/solid.blueprint" --thickness 4
```
will produce a partially smoothed 4-layer rendition of the solid as a
FTD blueprint.

![Smooth solid example](img/teapot%20smooth.jpg)

```
python convert.py "/path/to/some/solid.json" "/path/to/Constructables/some_donor.blueprint" "/path/to/Constructables/solid.blueprint" --thickness 4 --no-smoothing
```
will produce a coarse 4-layer rendition of the solid as a FTD blueprint.

![Coarse solid example](img/teapot%20coarse.jpg)

Tweaking
--------
Until there's a CLI option to do this (TODO), if you need to change the
way input coordinates are mapped to output coordinates (e.g. your model
is importing on its side or facing the wrong way), you can edit the part
of the script commented "transform the voxel coordinates".

Other limitations
-----------------
* Smoothing does not handle all possibilities (e.g. square corners are
  not handled), and is limited to 1-block slopes and corners.
* Converting a too large solid will take a long time and a lot of RAM.
* 32 layers maximum (color limit).

License
-------
Apache 2.0 (see LICENSE, or explainer [here](https://choosealicense.com/licenses/apache-2.0/)).