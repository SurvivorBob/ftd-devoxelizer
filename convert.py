#!/usr/bin/env python3
"""
Converts a voxel file into a FTD blueprint (see README).

License information:

   Copyright 2022 SurvivorBob <ftd-devoxelizer@survivorbob.xyz>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""

import json
import sys
import argparse

def neighbors(coords):
    return (
        (coords[0], coords[1], coords[2] - 1),
        (coords[0], coords[1], coords[2] + 1),
        (coords[0], coords[1] - 1, coords[2]),
        (coords[0], coords[1] + 1, coords[2]),
        (coords[0] - 1, coords[1], coords[2]),
        (coords[0] + 1, coords[1], coords[2]),
    )

single_block_id = 303
beam_ids = {
    2: 1050,
    3: 1051,
    4: 1052,
}
slope_id = 1055
corner_id = 1079
invcorner_id = 1060

def main():
    ap = argparse.ArgumentParser(description="Converts a JSON voxel encoding to a FTD blueprint, copying the author tag from a donor blueprint.")
    ap.add_argument("input_file", type=str, help="The input JSON voxel file (e.g. as produced by https://drububu.com/miscellaneous/voxelizer/?out=jso). Must be _solid_ (shell thickness 2048).")
    ap.add_argument("donor_blueprint", type=str, help="The donor blueprint from which to copy the author tag.")
    ap.add_argument("output_blueprint", type=str, help="The output file name for the blueprint to produce.")
    ap.add_argument("--thickness", type=int, default=2, help="Number of layers to retain (default 2, max 32).")
    ap.add_argument("--no-smoothing", action='store_true', help="Skips exterior smoothing if specified.")
    ap.add_argument("--no-beams", action='store_true', help="Skips beam consolidation if specified.")

    args = ap.parse_args(sys.argv[1:])

    if args.thickness < 1:
        print(f"Need to request at least one layer (requested {args.thickness})!")
        ap.print_usage()
        exit(-1)
    if args.thickness > 32:
        print(f"Tool supports only up to 32 layers thickness (requested {args.thickness})!")
        ap.print_usage()
        exit(-1)

    with open(args.input_file, mode="r") as input_file:
        voxel_data = json.load(input_file)

    with open(args.donor_blueprint, mode="r") as donor_blueprint_file:
        donor_blueprint = json.load(donor_blueprint_file)

    # preprocess the donor blueprint
    del donor_blueprint["Blueprint"]["VehicleData"]
    del donor_blueprint["Blueprint"]["CSI"]
    del donor_blueprint["SavedMaterialCost"]
    donor_blueprint["Blueprint"]["ContainedMaterialCost"] = 0.0
    donor_blueprint["ItemDictionary"]["303"] = "3cc75979-18ac-46c4-9a5b-25b327d99410" # single alloy block
    donor_blueprint["ItemDictionary"]["1050"] = "8f9dbf41-6c2d-4e7b-855d-b2432c6942a2" # 2 alloy block
    donor_blueprint["ItemDictionary"]["1051"] = "649f2aec-6f59-4157-ac01-0122ce2e6dad" # 3 alloy block
    donor_blueprint["ItemDictionary"]["1052"] = "9411e401-27da-4546-b805-3334f200f055" # 4 alloy block
    donor_blueprint["ItemDictionary"]["1055"] = "911fe222-f9b2-4892-9cd6-8b154d55b2aa" # 1 alloy slope
    donor_blueprint["ItemDictionary"]["1079"] = "a4b0d100-c480-4697-b606-489d80a6d376" # 1 alloy triangle corner
    donor_blueprint["ItemDictionary"]["1060"] = "95a626e6-f1b8-491a-aa31-8de5a2beb513" # 1 alloy triangle inverted
    donor_blueprint["Blueprint"]["BLP"] = []
    donor_blueprint["Blueprint"]["BLR"] = []
    donor_blueprint["Blueprint"]["BCI"] = []
    donor_blueprint["Blueprint"]["BlockIds"] = []
    donor_blueprint["Blueprint"]["TotalBlockCount"] = 0
    donor_blueprint["Blueprint"]["AliveCount"] = 0

    # +z == forward (rotation 0)
    # +y == up (rotation 10)
    # +x == right (rotation 1)
    x_min, x_max, y_min, y_max, z_min, z_max = 0, 0, 0, 0, 0, 0

    voxel_count = len(voxel_data["voxels"])

    src_width = int(voxel_data["dimension"][0]["width"]) + 1
    src_length = int(voxel_data["dimension"][0]["height"]) + 1
    src_height = int(voxel_data["dimension"][0]["depth"]) + 1

    src_x_center = int((src_width - 1) / 2)


    # transform the voxel coordinates
    voxels = set()
    for voxel in voxel_data["voxels"]:
        src_x, src_y, src_z = int(voxel["x"]), int(voxel["y"]), int(voxel["z"])
        dst_x = src_x - src_x_center
        dst_y = src_z
        dst_z = src_y

        x_min = min(dst_x, x_min)
        x_max = max(dst_x, x_max)
        y_min = min(dst_y, y_min)
        y_max = max(dst_y, y_max)
        z_min = min(dst_z, z_min)
        z_max = max(dst_z, z_max)

        voxels.add((dst_x, dst_y, dst_z))

    donor_blueprint["Blueprint"]["MaxCords"] = "{0},{1},{2}".format(x_max, y_max, z_max)
    donor_blueprint["Blueprint"]["MinCords"] = "{0},{1},{2}".format(x_min, y_min, z_min)

    # compute the depth of all the voxels (bellman-ford-like construction)
    # note: assumes a filled interior. TODO preprocessing required if given a shell
    voxel_depths = {}
    
    i = 0
    while True:
        no_change = True
        for v in voxels:
            minDist = sys.maxsize
            for n in neighbors(v):
                if n not in voxels:
                    minDist = 0
                elif n in voxel_depths:
                    minDist = min(voxel_depths[n], minDist)
            if minDist < sys.maxsize:
                if v not in voxel_depths:
                    no_change = False
                    voxel_depths[v] = minDist + 1
                elif voxel_depths[v] != minDist + 1:
                    no_change = False
                    voxel_depths[v] = minDist + 1
        if no_change:
            break
        
        print(f"voxel_depths {len(voxel_depths)} iter {i}")
        i += 1

    max_depth = max(voxel_depths.values())
    print(f"max depth: {max_depth}")

    # filter voxels by depth
    voxels_by_depth = {}
    for d in range(min(max_depth, args.thickness)):
        voxels_by_depth[d + 1] = set()
    
    for v, d in voxel_depths.items():
        if d in voxels_by_depth:
            voxels_by_depth[d].add(v)

    for d in voxels_by_depth.keys():
        print(f"depth {d}: {len(voxels_by_depth[d])}")

    def place_block(blueprint, v, bid, rot = 0, color = 0):
        blueprint["Blueprint"]["BLP"].append("{0},{1},{2}".format(v[0], v[1], v[2]))
        blueprint["Blueprint"]["BLR"].append(rot)
        blueprint["Blueprint"]["BCI"].append(color)
        blueprint["Blueprint"]["BlockIds"].append(bid)
        blueprint["Blueprint"]["TotalBlockCount"] += 1
        blueprint["Blueprint"]["AliveCount"] += 1

    # do corner erosion on outermost surface
    if not args.no_smoothing:
        outermost_voxels = voxels_by_depth[1]

        # bits:
        BITS_POSX = 0b00100000
        BITS_NEGX = 0b00010000
        BITS_POSY = 0b00001000
        BITS_NEGY = 0b00000100
        BITS_POSZ = 0b00000010
        BITS_NEGZ = 0b00000001

        def reverse_dirbits(bits):
            ret = 0
            if bits & BITS_POSX:
                ret = ret | BITS_NEGX
            if bits & BITS_NEGX:
                ret = ret | BITS_POSX
            if bits & BITS_POSY:
                ret = ret | BITS_NEGY
            if bits & BITS_NEGY:
                ret = ret | BITS_POSY
            if bits & BITS_POSZ:
                ret = ret | BITS_NEGZ
            if bits & BITS_NEGZ:
                ret = ret | BITS_POSZ
            return ret

        def neighbors_with_dirs(v):
            return (
                ((v[0] + 1, v[1], v[2]), BITS_POSX),
                ((v[0] - 1, v[1], v[2]), BITS_NEGX),
                ((v[0], v[1] + 1, v[2]), BITS_POSY),
                ((v[0], v[1] - 1, v[2]), BITS_NEGY),
                ((v[0], v[1], v[2] + 1), BITS_POSZ),
                ((v[0], v[1], v[2] - 1), BITS_NEGZ)
            )


        # returns a bitmap describing the orientation of neighboring empty space
        def empty_space_bitmap(v):
            if v not in voxels:
                return 0
            b = 0
            if (v[0] + 1, v[1], v[2]) not in voxels:
                b |= BITS_POSX
            if (v[0] - 1, v[1], v[2]) not in voxels:
                b |= BITS_NEGX
            if (v[0], v[1] + 1, v[2]) not in voxels:
                b |= BITS_POSY
            if (v[0], v[1] - 1, v[2]) not in voxels:
                b |= BITS_NEGY
            if (v[0], v[1], v[2] + 1) not in voxels:
                b |= BITS_POSZ
            if (v[0], v[1], v[2] - 1) not in voxels:
                b |= BITS_NEGZ
            return b

        def adjacency_bitmap(v, other_set):
            if v not in voxels:
                return 0
            b = 0
            if (v[0] + 1, v[1], v[2]) in other_set:
                b |= BITS_POSX
            if (v[0] - 1, v[1], v[2]) in other_set:
                b |= BITS_NEGX
            if (v[0], v[1] + 1, v[2]) in other_set:
                b |= BITS_POSY
            if (v[0], v[1] - 1, v[2]) in other_set:
                b |= BITS_NEGY
            if (v[0], v[1], v[2] + 1) in other_set:
                b |= BITS_POSZ
            if (v[0], v[1], v[2] - 1) in other_set:
                b |= BITS_NEGZ
            return b
        
        # returns a (block type, orientation, check for layer 2, implies square face in direction) based on empty space
        # this is a table lookup because deriving this computationally sounds horrific
        cornerness_map = {
            # slopes
            BITS_POSX | BITS_POSY: (slope_id, 1, BITS_POSZ | BITS_NEGZ, BITS_NEGX | BITS_NEGY),
            BITS_POSX | BITS_NEGY: (slope_id, 13, BITS_POSZ | BITS_NEGZ, BITS_NEGX | BITS_POSY),
            BITS_POSX | BITS_POSZ: (slope_id, 16, BITS_POSY | BITS_NEGY, BITS_NEGX | BITS_NEGZ),
            BITS_POSX | BITS_NEGZ: (slope_id, 22, BITS_POSY | BITS_NEGY, BITS_NEGX | BITS_POSZ),
            BITS_NEGX | BITS_POSY: (slope_id, 3, BITS_POSZ | BITS_NEGZ, BITS_POSX | BITS_NEGY),
            BITS_NEGX | BITS_NEGY: (slope_id, 15, BITS_POSZ | BITS_NEGZ, BITS_POSX | BITS_POSY),
            BITS_NEGX | BITS_POSZ: (slope_id, 21, BITS_POSY | BITS_NEGY, BITS_POSX | BITS_NEGZ),
            BITS_NEGX | BITS_NEGZ: (slope_id, 19, BITS_POSY | BITS_NEGY, BITS_POSX | BITS_POSZ),
            BITS_POSY | BITS_POSZ: (slope_id, 0, BITS_POSX | BITS_NEGX, BITS_NEGY | BITS_NEGZ),
            BITS_POSY | BITS_NEGZ: (slope_id, 2, BITS_POSX | BITS_NEGX, BITS_NEGY | BITS_POSZ),
            BITS_NEGY | BITS_POSZ: (slope_id, 12, BITS_POSX | BITS_NEGX, BITS_POSY | BITS_NEGZ),
            BITS_NEGY | BITS_NEGZ: (slope_id, 14, BITS_POSX | BITS_NEGX, BITS_POSY | BITS_POSZ),

            # corners
            BITS_POSX | BITS_POSY | BITS_POSZ: (corner_id, 0, BITS_NEGX | BITS_NEGY | BITS_NEGZ, 0),
            BITS_POSX | BITS_POSY | BITS_NEGZ: (corner_id, 1, BITS_NEGX | BITS_NEGY | BITS_POSZ, 0),
            BITS_POSX | BITS_NEGY | BITS_POSZ: (corner_id, 16, BITS_NEGX | BITS_POSY | BITS_NEGZ, 0),
            BITS_POSX | BITS_NEGY | BITS_NEGZ: (corner_id, 22, BITS_NEGX | BITS_POSY | BITS_POSZ, 0),
            BITS_NEGX | BITS_POSY | BITS_POSZ: (corner_id, 3, BITS_POSX | BITS_NEGY | BITS_NEGZ, 0),
            BITS_NEGX | BITS_POSY | BITS_NEGZ: (corner_id, 2, BITS_POSX | BITS_NEGY | BITS_POSZ, 0),
            BITS_NEGX | BITS_NEGY | BITS_POSZ: (corner_id, 21, BITS_POSX | BITS_POSY | BITS_NEGZ, 0),
            BITS_NEGX | BITS_NEGY | BITS_NEGZ: (corner_id, 19, BITS_POSX | BITS_POSY | BITS_POSZ, 0),

        }
        def cornerness(v):
            return cornerness_map.get(empty_space_bitmap(v), (None, None, None, None))

        # maps (empty space + slopes + corners bitmap to inverse corners)
        inverse_cornerness_map = {
            BITS_POSX | BITS_POSY | BITS_POSZ: (invcorner_id, 0, 0, BITS_NEGX | BITS_NEGY | BITS_NEGZ),
            BITS_POSX | BITS_POSY | BITS_NEGZ: (invcorner_id, 1, 0, BITS_NEGX | BITS_NEGY | BITS_POSZ),
            BITS_POSX | BITS_NEGY | BITS_POSZ: (invcorner_id, 16, 0, BITS_NEGX | BITS_POSY | BITS_NEGZ),
            BITS_POSX | BITS_NEGY | BITS_NEGZ: (invcorner_id, 22, 0, BITS_NEGX | BITS_POSY | BITS_POSZ),
            BITS_NEGX | BITS_POSY | BITS_POSZ: (invcorner_id, 3, 0, BITS_POSX | BITS_NEGY | BITS_NEGZ),
            BITS_NEGX | BITS_POSY | BITS_NEGZ: (invcorner_id, 2, 0, BITS_POSX | BITS_NEGY | BITS_POSZ),
            BITS_NEGX | BITS_NEGY | BITS_POSZ: (invcorner_id, 21, 0, BITS_POSX | BITS_POSY | BITS_NEGZ),
            BITS_NEGX | BITS_NEGY | BITS_NEGZ: (invcorner_id, 19, 0, BITS_POSX | BITS_POSY | BITS_POSZ),

            BITS_POSX | BITS_POSY | BITS_POSZ | BITS_NEGZ: (slope_id, 1, 0, BITS_NEGX | BITS_NEGY),
            BITS_POSX | BITS_NEGY | BITS_POSZ | BITS_NEGZ: (slope_id, 13, 0, BITS_NEGX | BITS_POSY),
            BITS_POSX | BITS_POSZ | BITS_POSY | BITS_NEGY: (slope_id, 16, 0, BITS_NEGX | BITS_NEGZ),
            BITS_POSX | BITS_NEGZ | BITS_POSY | BITS_NEGY: (slope_id, 22, 0, BITS_NEGX | BITS_POSZ),
            BITS_NEGX | BITS_POSY | BITS_POSZ | BITS_NEGZ: (slope_id, 3, 0, BITS_POSX | BITS_NEGY),
            BITS_NEGX | BITS_NEGY | BITS_POSZ | BITS_NEGZ: (slope_id, 15, 0, BITS_POSX | BITS_POSY),
            BITS_NEGX | BITS_POSZ | BITS_POSY | BITS_NEGY: (slope_id, 21, 0, BITS_POSX | BITS_NEGZ),
            BITS_NEGX | BITS_NEGZ | BITS_POSY | BITS_NEGY: (slope_id, 19, 0, BITS_POSX | BITS_POSZ),
            BITS_POSY | BITS_POSZ | BITS_POSX | BITS_NEGX: (slope_id, 0, 0, BITS_NEGY | BITS_NEGZ),
            BITS_POSY | BITS_NEGZ | BITS_POSX | BITS_NEGX: (slope_id, 2, 0, BITS_NEGY | BITS_POSZ),
            BITS_NEGY | BITS_POSZ | BITS_POSX | BITS_NEGX: (slope_id, 12, 0, BITS_POSY | BITS_NEGZ),
            BITS_NEGY | BITS_NEGZ | BITS_POSX | BITS_NEGX: (slope_id, 14, 0, BITS_POSY | BITS_POSZ),
        }
        def invcornerness(v, corners, slopes, square_in_dir_map):
            corner_adjacency = adjacency_bitmap(v, corners) | adjacency_bitmap(v, slopes)

            for n, dir_to_n in neighbors_with_dirs(v):
                dir_to_v = reverse_dirbits(dir_to_n)
                # if v needs to be square in the direction from n, then we can't consider n as a slope or corner
                if dir_to_v & square_in_dir_map.get(n, 0):
                    corner_adjacency = corner_adjacency & (0xFF ^ dir_to_n)

            return inverse_cornerness_map.get(empty_space_bitmap(v) | corner_adjacency, (None, None, None, None)) + (corner_adjacency,)

        placed_corners = set()
        placed_slopes = set()

        square_in_dir_for_voxels = {}

        placed_invcorners = set()

        for v in outermost_voxels:
            block_id, rot, to_test, square_in_dir = cornerness(v)
            placing_block = False
            if block_id is not None and rot is not None:
                if block_id == slope_id and adjacency_bitmap(v, voxels_by_depth[2]) & to_test == 0:
                    placed_slopes.add(v)
                    placing_block = True
                elif block_id == corner_id and adjacency_bitmap(v, voxels_by_depth[2]) & to_test == 0:
                    placed_corners.add(v)
                    placing_block = True

                if placing_block:
                    square_in_dir_for_voxels[v] = square_in_dir_for_voxels.get(v, 0) | square_in_dir
                    place_block(donor_blueprint, v, block_id, rot, 0)
        
        voxels_by_depth[1] = outermost_voxels - placed_corners - placed_slopes
        outermost_voxels = voxels_by_depth[1]

        for v in outermost_voxels:
            block_id, rot, to_test, square_in_dir, corner_adjacency = invcornerness(v, placed_corners, placed_slopes, square_in_dir_for_voxels)
            placing_block = False
            if block_id is not None and rot is not None:
                if block_id == invcorner_id and corner_adjacency > 0 and adjacency_bitmap(v, voxels_by_depth[2]) & to_test == 0:
                    placed_invcorners.add(v)
                    placing_block = True
                elif block_id == slope_id and adjacency_bitmap(v, voxels_by_depth[2]) & to_test == 0:
                    placed_slopes.add(v)
                    placing_block = True

                if placing_block:
                    square_in_dir_for_voxels[v] = square_in_dir_for_voxels.get(v, 0) | square_in_dir
                    place_block(donor_blueprint, v, block_id, rot, 0)
        
        voxels_by_depth[1] = outermost_voxels - placed_invcorners - placed_slopes
        outermost_voxels = voxels_by_depth[1]

        print(f"d 1: created {len(placed_slopes)} slopes, {len(placed_corners)} triangle corners, {len(placed_invcorners)} invcorners")
        print(f"d 1: {len(voxels_by_depth[1])} voxels left")

    # consolidate beams in z and y axes
    if not args.no_beams:

        # consolidate beams (z axis, rotation 0)
        for d in voxels_by_depth.keys():
            visited_voxels = set()
            n_beams = 0
            for v in voxels_by_depth[d]:
                beam_fd = 0
                beam_bk = 0
                if v not in visited_voxels:
                    # scan forward
                    b = v
                    while beam_fd + 1 < 4:
                        b = (b[0], b[1], b[2] + 1)
                        if b not in visited_voxels and b in voxels_by_depth[d]:
                            beam_fd += 1
                        else:
                            break
                    # scan backward
                    b = v
                    while beam_bk + 1 < 4 - beam_fd:
                        b = (b[0], b[1], b[2] - 1)
                        if b not in visited_voxels and b in voxels_by_depth[d]:
                            beam_bk += 1
                        else:
                            break

                    # find real origin point of beam
                    beam_origin = (v[0], v[1], v[2] - beam_bk)
                    beam_length = 1 + beam_fd + beam_bk

                    if beam_length > 1:
                        block_id = beam_ids[beam_length]
                        place_block(donor_blueprint, beam_origin, block_id, 0, d - 1)
                        n_beams += 1
                        for i in range(beam_length):
                            pt = (beam_origin[0], beam_origin[1], beam_origin[2] + i)
                            if pt in visited_voxels:
                                print(f"d {d}: pt {pt} already visited???")
                            visited_voxels.add(pt)
            print(f"d {d}: turned {len(visited_voxels)} voxels into {n_beams} beams")
            voxels_by_depth[d] = voxels_by_depth[d] - visited_voxels
            print(f"d {d}: {len(voxels_by_depth[d])} voxels left")

        # consolidate beams (y axis, rotation 10)

        for d in voxels_by_depth.keys():
            visited_voxels = set()
            n_beams = 0
            for v in voxels_by_depth[d]:
                beam_fd = 0
                beam_bk = 0
                if v not in visited_voxels:
                    # scan forward
                    b = v
                    while beam_fd + 1 < 4:
                        b = (b[0], b[1] + 1, b[2])
                        if b not in visited_voxels and b in voxels_by_depth[d]:
                            beam_fd += 1
                        else:
                            break
                    # scan backward
                    b = v
                    while beam_bk + 1 < 4 - beam_fd:
                        b = (b[0], b[1] - 1, b[2])
                        if b not in visited_voxels and b in voxels_by_depth[d]:
                            beam_bk += 1
                        else:
                            break

                    # find real origin point of beam
                    beam_origin = (v[0], v[1] - beam_bk, v[2])
                    beam_length = 1 + beam_fd + beam_bk

                    if beam_length > 1:
                        block_id = beam_ids[beam_length]
                        place_block(donor_blueprint, beam_origin, block_id, 10, d - 1)
                        n_beams += 1
                        for i in range(beam_length):
                            pt = (beam_origin[0], beam_origin[1] + i, beam_origin[2])
                            if pt in visited_voxels:
                                print(f"d {d}: pt {pt} already visited???")
                            visited_voxels.add(pt)
            print(f"d {d}: turned {len(visited_voxels)} voxels into {n_beams} beams")
            voxels_by_depth[d] = voxels_by_depth[d] - visited_voxels
            print(f"d {d}: {len(voxels_by_depth[d])} voxels left")

        # TODO: consolidate x beams (probably need to cut model in half lol)

    # place remaining singleton blocks
    for d in voxels_by_depth.keys():
        for v in voxels_by_depth[d]:
            place_block(donor_blueprint, v, single_block_id, 0, d - 1)

    # do final fixup here

    donor_blueprint["Blueprint"]["BlockState"] = "=0,{0}".format(donor_blueprint["Blueprint"]["TotalBlockCount"])

    donor_blueprint["SavedTotalBlockCount"] = donor_blueprint["Blueprint"]["TotalBlockCount"]
    # donor_blueprint["SavedMaterialCost"] = voxel_count * 5.0
    print("saving...")
    with open(args.output_blueprint, mode="w") as output_blueprint_file:
        json.dump(donor_blueprint, output_blueprint_file)

    print("all done!")

if __name__ == "__main__":
    main()