#!/usr/bin/env python3

import argparse

COV_COL_IDX = 5
DEPTH_COL_IDX = 6

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find the accession number associated with the highest bbmap median_fold for each virus species."
    )
    parser.add_argument(
        "-bbmap_covstats", required=True, type=str, help="bbmap covstats output file"
    )
    parser.add_argument("-b", required=True, type=str, help="sample basename")
    parser.add_argument(
        "-m",
        type=int,
        default=3,
        help="Minimum mean depth for a reference to be selected (Default 3)",
    )
    parser.add_argument(
        "-p",
        type=int,
        default=70,
        help="Minimum % coverage for a reference to be selected (Default 70)",
    )
    args = parser.parse_args()

    init_ref_candidates = []
    to_sort = []
    inf = open(args.bbmap_covstats, "r").readlines()
    header = inf[0]
    rec = inf[1:]

    # add to list any reference that passes the median coverage and mininum covered percent threshold
    for i in rec:
        if len(i) > 0:
            temp = i.split("\t")
            if (
                float(temp[DEPTH_COL_IDX]) >= args.m
                and float(temp[COV_COL_IDX]) >= args.p
            ):
                print(
                    "Found reference above mean depth (%s, %s) and mean coverage (%s, %s)"
                    % (
                        args.m,
                        float(temp[DEPTH_COL_IDX]),
                        args.p,
                        float(temp[COV_COL_IDX]),
                    )
                )
                # init_ref_candidates.append(temp[0])
                to_sort.append(temp)

    sorted_by_coverage_depth = sorted(to_sort, key=lambda x: (float(x[COV_COL_IDX]), float(x[DEPTH_COL_IDX])), reverse=True)
    init_ref_candidates = [entry[0] for entry in sorted_by_coverage_depth]

    # create a dictionary of references to be used for consensus calling
    init_ref_header = {}
    for i in init_ref_candidates:
        # add .split(' ')[0] to get just the accession
        acc = i.split(" ")[0]
        # get header tag
        tag = i.split(" ")[1]
        # get header info
        info = " ".join(i.split(" ")[2:])
        if not tag in init_ref_header:
            init_ref_header[tag] = [acc, info]

    # if reference(s) selected, output relevant info to file
    if init_ref_header:
        output_file_name = args.b + "_refs.tsv"
        output_file = open(output_file_name, "a+")
        for i in init_ref_header:
            # output format: reference accession <tab> reference header tag <tab> reference header info
            output_file.write(
                str(init_ref_header[i][0])
                + "\t"
                + str(i)
                + "\t"
                + str(init_ref_header[i][1])
                + "\n"
            )
        output_file.close()

    # no reference selected, output reference with highest covered percent and read distribution info
    else:
        output_file_name = args.b + "_failed_assembly.tsv"
        output_file = open(output_file_name, "w")

        output_text = str(
            sorted(rec, key=lambda line: float(line.split("\t")[4]), reverse=True)[0]
        )
        output_file.write(header)
        output_file.write("\n")
        output_file.write(output_text)
        output_file.write("\n")
        output_file.close()
