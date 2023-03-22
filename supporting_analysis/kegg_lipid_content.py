""" A script that looks into the lipid represntation in KEGG using the lipids BRITE json file downloaded in 02/14/2023"""

import json
from pathlib import Path
import os


def count_children_in_tree(lipid_dict, total_kegg_lipids):
    if "children" in lipid_dict:
        next_dict=lipid_dict["children"]
        for level_kid in next_dict:
            count_children_in_tree(level_kid, total_kegg_lipids)
    else:
        total_kegg_lipids.add(lipid_dict["name"])
    return total_kegg_lipids


def main():
    """ The main function"""
    util_path=Path("./util_files")
    kegg_lipid_brite_dict=json.load(open(os.path.join(util_path, "br08002.json")))
    total_kegg_lipids=set()
    total_kegg_lipids=count_children_in_tree(kegg_lipid_brite_dict, total_kegg_lipids)

    print(len(total_kegg_lipids))




if __name__ == "__main__":

    main()

    print("Done")
