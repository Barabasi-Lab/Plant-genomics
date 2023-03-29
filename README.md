
![Logo](https://github.com/Barabasi-Lab/Plant-genomics/blob/main/images/NetSci_Logo.png)



# Genomics-based annotations help unveil the molecular composition of edible plants

Nutrition and wellbeing take a central role in today’s high pace life, but how much do we really know about the food we eat? Here, we harness existing metabolic knowledge encrypted in staple food ingredients’ genome to help us explore the composition of raw edible plants. We first show the benefit and value of looking into genome-associated functional annotations on a wide scale. Next, we rely on new experimental data to develop a framework that helps us reveal new, potentially bioactive compounds in staple food ingredients. This has significance in, first, extending current food composition knowledge and second in discovering newly detected bioactive compounds, shedding light on the potential impacts of common food ingredients beyond their nutritional value as described in food labels. Finally, we show that staple foods that are already included in our daily diets might have the potential to be ‘superfoods’ that can contribute to our wellbeing.

Here we explore to what degree, existing genome-associated metabolic annotations can offer a valuable resource to deepen the knowledge of food composition. To do so and focus on the edible parts of plants, we rely on metabologenomics, integrating genomics and metabolomics, used in the past to discover novel natural products. 

## Installation

To use the codein this repository please install the requirments in ```requirments.txt```



Use pip install -r requirements.txt to install the related packages.

## Code and Data

### Code

``` plant_analysis```

**General plant analysis**

PathEnrichmentPerOrg_bonferroni.py:

**Corn specific analysis**

corn_analysis.py:


structure_analysisCorn.py:

Calculates the similarity between corn related compounds and creates the clustermap presented in figure 5.

``` thermodynamic_feasibility_analysis ```

feasibleReactionsPerOrg.py:


kineticsPerformance_fb.py:


``` supporting_analysis ```

coconutdb_comparison.py:

kegg_lipid_content.py:

pca_approch.ipynb:

plant_family.py:

### Data

All the data files are ocated in the data folder.

Metabolomics data (Metabolon_Data.xlsx) - Contains the metabolomics experiments results for the plants presented in the paper.

Metabolic annotations master table  (plant_masterTable_100621_ds1.csv) - Contains all the annotations collected for the plants presented in the paper.

Table_s2- A comparison of first block InchIKeys from coconutDB, our predicted to accumulate compounds, and experimentally detected compounds.

Thermodynamic feasibility score table (kinetics_all_plants_df_fb_avg.csv) - Atable summarizing the thermodynamic feasibility approach results including the score, is it found in the experiments, first block inchi key, the group of full inchi keys for that first block representation, compound mass and the number of reactions the compound is in.
