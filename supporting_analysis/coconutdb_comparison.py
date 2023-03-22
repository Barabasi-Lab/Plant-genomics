""" This script compares the first block inchis for the metabolomics experimnets"""
import matplotlib
matplotlib.use('tkagg')
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
from pathlib import Path
import pickle
import pubchempy as pcp



def parse_plant_master_table(master_table_file, fb2smiles):

    master_table_df=pd.read_csv(master_table_file)
    for col in master_table_df:
        col_list=master_table_df[col].tolist()
        fb=col_list[1]
        if "fbloack" in fb:
            continue
        smiles=col_list[2]
        if fb not in fb2smiles:
            fb2smiles[fb]=smiles
    return fb2smiles

def compile_coconutdb_first_block_set(ccdb_file_path, fb2smiles):
    """Loads the coconutDB file used for comparison. The database dump is current to the 2022 release."""
    ccdb_df=pd.read_csv(ccdb_file_path)
    # only use the recent release
    ccdb_df_recent=ccdb_df[ccdb_df["year"]!="2020"]
    for i,row in ccdb_df_recent.iterrows():
        fb=row["FB"]
        smiles=row["clean_smiles"]
        if fb not in fb2smiles:
            fb2smiles[fb]=smiles
    ccdb_fb_set=set(ccdb_df_recent["FB"].unique().tolist())
    return ccdb_fb_set, fb2smiles

def convert2fb(iklist):
    """ converts a full key list to first block list"""
    fb2ik={}
    fblist=set()
    for ik in iklist:
        ikl=ik.split('-')
        if ikl[0].startswith('XXFWWUCLFPRWQI'):
            pass
        fblist.add(ikl[0])
        fb2ik[ikl[0]]=ik
    return fblist, fb2ik

def compile_experimental_first_block_set(data_file):
    """Loades the our plant data master table and gets the complete list of
    compounds experimentally detected.
    """


    # work on the plants we have experiments for
    orgs=pickle.load(open(data_file,'rb'))
    totonlyexpri=set()
    plantlist=set()
    totexp=set()
    # get big M for the overlaps with experiments
    for org in orgs:
        if 'experiment' in orgs[org]:
            plantlist.add(org)
            for atr in orgs[org]:
                if atr in orgs[org]:
                    totonlyexpri = totonlyexpri | orgs[org][atr]
                    if 'experiment' in atr:
                        totexp=totexp|orgs[org][atr]
    experiment_fb_set, fb2ik=convert2fb(totexp)
    return experiment_fb_set



def compile_the_predicted_compound_first_block_set(thermodynamics_results_file_path):
    """ Parses and compiles the compounds predicted to accumulate by our thermodynamic feaasibility approach

    """
    thermodynamics_predictions_df=pd.read_csv(thermodynamics_results_file_path)
    thdb_predicted_to_accumulate_df=thermodynamics_predictions_df[thermodynamics_predictions_df["sum"]>0]

    return set(thdb_predicted_to_accumulate_df["inchi key"].unique().tolist())

def plot_3_way_venn(set1,set2,set3, label1, label2, label3, filename):
    """Plots a three way venn for given sets"""
    plt.figure(figsize = (10, 10))
    venn3([set1, set2, set3], (label1,label2,label3),
        set_colors=('lightskyblue', 'darkorchid', 'goldenrod'),
        alpha=0.8)
    plt.savefig(filename, dpi=600, bbox_inches = 'tight')

def compile_result_table(coconutdb_fb_set, predicted_fb_set, experiments_fb_set, fb2smiles):
    """Compiles a results table for this analysis"""
    pooled_fb_set=coconutdb_fb_set|predicted_fb_set|experiments_fb_set
    fb_wo_smiles=set()
    results_dict={}
    for fb in pooled_fb_set:
        results_dict[fb]={}
        if fb in fb2smiles:
            results_dict[fb]["smiles"]=fb2smiles[fb]
        else:
            fb_wo_smiles.add(fb)
        if fb in coconutdb_fb_set:
            results_dict[fb]["CoconutDB"]=1
        else:
            results_dict[fb]["CoconutDB"]=0
        if fb in predicted_fb_set:
            results_dict[fb]["Predicted to accumulate"]=1
        else:
            results_dict[fb]["Predicted to accumulate"]=0
        if fb in experiments_fb_set:
            results_dict[fb]["Experimentally detected"]=1
        else:
            results_dict[fb]["Experimentally detected"]=0

    results_df=pd.DataFrame.from_dict(results_dict)
    results_df=results_df.transpose()

    return results_df

def main():
    """The main function"""
    fb2smiles={}
    util_path=Path("./util_files")
    results_path=Path("./results/")

    coconutdb_fb_set, fb2smiles=compile_coconutdb_first_block_set(os.path.join(util_path,"Coconut-NDM_Overlap.csv"), fb2smiles)
    predicted_fb_set=compile_the_predicted_compound_first_block_set(os.path.join(util_path,"kinetics_all_Plants_df_fb_avg.csv"))
    experiments_fb_set=compile_experimental_first_block_set(os.path.join(util_path,"orgkeys_new.pkl"))
    fb2smiles=parse_plant_master_table(os.path.join(util_path,"plant_masterTable_100621.csv"), fb2smiles)
    plot_3_way_venn(coconutdb_fb_set,predicted_fb_set,experiments_fb_set,"CoconutDB","Predicted","Experimentally detected",os.path.join(results_path,"coconutdb_comparison_venn.svg"))

    experimental_not_in_coconut_and_not_in_pred=experiments_fb_set-coconutdb_fb_set-predicted_fb_set
    print(len(experimental_not_in_coconut_and_not_in_pred))
    results_df=compile_result_table(coconutdb_fb_set, predicted_fb_set, experiments_fb_set, fb2smiles)
    results_df.to_csv(os.path.join(results_path,"table_s2.csv"))

if __name__ == "__main__":

    main()

    print("Done")
