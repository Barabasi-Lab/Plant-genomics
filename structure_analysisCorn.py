''' This script performs a structural similarity and distance analysis on the list of smiles.'''



# imports
import sys

import pandas as pd
import numpy as np
import pickle
import time

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

#from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw

from functions_structural_similarity import get_binary_representation, get_count_representation, get_bit_info

IPythonConsole.ipython_useSVG=True
#from rdkit import DataStructs

from pandarallel import pandarallel
from tqdm import tqdm

from mol2vec.features import mol2alt_sentence, MolSentence, DfVec, sentences2vec, mol2sentence
from mol2vec.helpers import depict_identifier, mol_to_svg, IdentifierTable, plot_2D_vectors
from scipy.stats import hypergeom
import collections
import re
pandarallel.initialize(progress_bar = True,nb_workers=16)
tqdm.pandas()
import json
import umap
import umap.plot
from scipy.spatial.distance import pdist, squareform
from statsmodels.distributions.empirical_distribution import ECDF
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import joblib
from IPython.core.display import HTML
from IPython.display import SVG
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
PandasTools.RenderImagesInAllDataFrames(images=True)
from matplotlib_venn import venn2
import matplotlib
import itertools
from itertools import compress
from statsmodels.stats.multitest import fdrcorrection
import pubchempy as pcp
from bioservices import *





# ----------------- functions -----------------
def cleanhtml(raw_html):
    # html tag cleaner copied from stack overflow, cleans non tag entities as well
  cleanr = re.compile('<.*?>|&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-f]{1,6});')
  cleantext = re.sub(cleanr, '', raw_html)
  return cleantext
def getName(inchikey):
    ''' Gets the name and isomeric smiles string from pubchem'''
    while True:
        try:
            revcomp = pcp.get_compounds(inchikey, 'inchikey')
            if revcomp:
                synonyms=revcomp[0].synonyms
                #name=revcomp[0].iupac_name
                if synonyms:
                    name=synonyms[0]
                else:
                    name = revcomp[0].iupac_name
                # update mdm
                mdm[inchikey]['name']=name

                break
            else:
                # look in chebi
                if 'chebi' in mdm[inchikey]:
                    ccomp = ch.getLiteEntity(mdm[inchikey]['chebi'][0])
                    chcomp=ch.getCompleteEntity(mdm[inchikey]['chebi'][0])
                    if chcomp:
                        smiles=chcomp.smiles
                        name=ccomp[0].chebiAsciiName
                        mdm[inchikey]['name'] = name

                        break
        except Exception as ex:
            # If an error is raised it will print it while the code runs.
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print(exc_type, exc_tb.tb_lineno)
            print('Retrying in 0.5 seconds. Exception in connecting')
            time.sleep(0.5)
    return name

def getSmiles(inchikey):
    ''' Gets the name and isomeric smiles string from pubchem'''
    while True:
        try:
            revcomp = pcp.get_compounds(inchikey, 'inchikey')
            if revcomp:
                smiles=revcomp[0].isomeric_smiles

                # update mdm

                mdm[inchikey]['smiles']=smiles
                break
            else:
                # look in chebi
                if 'chebi' in mdm[inchikey]:
                    ccomp = ch.getLiteEntity(mdm[inchikey]['chebi'][0])
                    chcomp=ch.getCompleteEntity(mdm[inchikey]['chebi'][0])
                    if chcomp:
                        smiles=chcomp.smiles


                        mdm[inchikey]['smiles'] = smiles
                        break
        except Exception as ex:
            # If an error is raised it will print it while the code runs.
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print(exc_type, exc_tb.tb_lineno)
            print('Retrying in 0.5 seconds. Exception in connecting')
            time.sleep(0.5)
    return smiles
def addMolecule(dataframe):
    ''' Adds a rdkit molecule represntation to a datframe'''
    PandasTools.AddMoleculeColumnToFrame(dataframe, 'smiles', 'Molecule', includeFingerprints=True)


def getFpBitVec(dataframe):
    ''' generates a 8,192 bit vector for a radius 3 morgan fingerprint'''

    dataframe['binary_f'] = dataframe['smiles'].progress_apply(get_binary_representation, radius=3,
                                                                                   nBits=8192)
    dataframe['count_f'] = dataframe['smiles'].progress_apply(get_count_representation, radius=3,
                                                                                  nBits=8192)
    dataframe['len_bit'] = dataframe['smiles'].progress_apply(get_binary_representation, radius=3,
                                                                                  nBits=8192).apply(np.sum)
    dataframe['bit_info'] = dataframe['smiles'].progress_apply(get_bit_info, radius=3, nBits=8192)


def getEmbeddings(dataframe):
    ''' create the embeddings for a dataframe bit representation'''
    Embeddings = np.array([np.array(xi) for xi in dataframe['binary_f']])
    nc = Embeddings.shape[0]
    nf = Embeddings.shape[1]
    kf = Embeddings.sum(axis=0)
    return Embeddings


def createUMAP(embeds):
    reducer = umap.UMAP(n_components=2, metric='jaccard')
    umapcoord = reducer.fit(embeds)
    return umapcoord


def getJacSimMtx(embeds, dataframe):
    matdist = squareform(pdist(embeds, 'jaccard'))
    matsim = 1 - matdist
    dfsim = pd.DataFrame(data=matsim, index=dataframe['name'], columns=dataframe['name'])
    return matsim,matdist


# ---------------- Main -----------------------
ch=ChEBI()

''' I will have three scenarios:
1. General - all molecules that were classified to all classes from all organisms
2. Plant level = all the molecules classified  to all classes in the plant
3. plant_class level - molecules per class per plant

How will this be visualized?'''

# part 1 - create the input for the structural analysis
orgs=pickle.load(open('orgs_classified_orgkeys.pkl','rb'))

# I need the smiles - I can find them in ndm and mdm
# read mdm
print('loading MDM')
'''mdm={}
with open('MDM_021721.json') as mdmin:
    mdm=json.load(mdmin)
'''
ndmdf=pd.read_csv('NDM_MasterJan21.csv')
ndm=ndmdf.set_index('InChIKey').to_dict()
mdm=pickle.load(open('mdm_040621.pkl','rb'))
ik2type=pickle.load(open('ik2prim_sec.pkl','rb'))
df=pd.DataFrame()
cp1=['dodgerblue','red','grey']
# create the dataframe for analysis
totkeys=set()
totsmiles=set()
ik2name=pickle.load(open('ik2iupacname.pkl','rb'))
ik2isosmiles=pickle.load(open('ik2isosmiles.pkl','rb'))
'''for org in tqdm(orgs):
    for cls in tqdm(orgs[org]):
        for ik in orgs[org][cls]:
            totkeys.add(ik)
            if ik in mdm:
                if ik in ik2name and ik in ik2isosmiles:
                    name=ik2name[ik]
                    smiles=ik2isosmiles[ik]
                else:
                    # need to get a name
                    smiles,name=getInfo(ik)
                    totsmiles.add(smiles)
                    ik2name[ik]=name
                    ik2isosmiles[ik]=smiles


            values_to_add = {'Class': cls, 'Inchi key': ik, 'smiles': smiles, 'Name': name, 'Plant': org}
            row_to_add = pd.Series(values_to_add, name=org)
            df = df.append(row_to_add)'''
# save the new mdmd
'''pickle.dump( mdm, open( "mdm_040621.pkl", "wb" ) )
pickle.dump( ik2name, open( "ik2iupacname.pkl", "wb" ) )
pickle.dump( ik2isosmiles, open( "ik2isosmiles.pkl", "wb" ) )'''
df=pd.read_csv('exp_gen_corn.csv')
# the input is ready - start analysis on the general list of keys
addMolecule(df)
# remove all the lines that do not have a molecule
# split the dataframe into smaller chunks
dfparts=np.array_split(df, 20)
new_df=pd.DataFrame()
getFpBitVec(df)

genembeds=getEmbeddings(df)
map=createUMAP(genembeds)
matsim,matdist=getJacSimMtx(genembeds,df)

#visualize everything

# count according to classes
classes=df['type'].unique()
ordered_classes=sorted(classes)
counterST=collections.Counter(df['type'])
counterSTkey=sorted(counterST, key=counterST.get, reverse=True)
counttypes=[counterST[k] for k in counterSTkey]
print([(k, counterST[k]) for k in counterSTkey])


'''plt.figure(figsize = (20, 20))
p=umap.plot.points(map)
plt.tight_layout()
plt.show()'''



plt.figure(figsize = (10, 10))
sns.set(font_scale=2)
#sns.color_palette("viridis", as_cmap=True)
sns.set_palette(sns.color_palette(cp1))
#sns.color_palette('viridis', n_colors=len(ordered_classes))
ax = sns.barplot(x = counttypes, y=counterSTkey, palette=sns.color_palette(cp1), edgecolor='gray')
plt.xlabel('Inchi key count')
plt.ylabel('Metabolism')
plt.tight_layout()
plt.show()


lut = dict(zip(set(counterSTkey), sns.color_palette(cp1)))
row_colors = list(df['type'].map(lut))


# plot the jaccard similarity
'''plt.figure(figsize = (20, 10))
sns.color_palette('viridis', n_colors=len(ordered_classes))
for lb in counterSTkey:
    msel=squareform(matdist[df['Class']==lb,:][:,df['Class']==lb])
    sns.distplot(1-msel, label=lb, hist=False,color=lut[lb])

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=5, mode="expand", borderaxespad=0., fontsize=12,markerscale=2.)
plt.xlabel('Jaccard Similarity')
plt.tight_layout()
plt.show()

plt.figure(figsize=(20, 10))

sns.color_palette('viridis', n_colors=len(counterSTkey))
for target in counterSTkey:
    sns.distplot(df.loc[df['Class'] == target, 'len_bit'], label=target, color=lut[target])

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=5, mode="expand", borderaxespad=0., fontsize=12,markerscale=2.)
plt.xlabel('Active Bits')
plt.ylabel('PDF')
plt.title('Morgan Fingerprint 8192 bits, Radius 3')
plt.tight_layout()
plt.show()'''

# boxplot
similarityvec = []
outcomvec = []

for lb in counterSTkey:
    msel = squareform(matdist[df['type'] == lb, :][:, df['type'] == lb])
    similarityvec.append(1 - msel)
    outcomvec.append([lb for ind in range(len(msel))])

similarityvec = list(itertools.chain(*similarityvec))
outcomvec = list(itertools.chain(*outcomvec))

plotdf = pd.DataFrame.from_dict({'Jaccard': similarityvec, 'Outcome': outcomvec})

plt.figure(figsize=(10, 8))
sns.set(font_scale=2)
sns.set_palette(sns.color_palette(cp1))
ax = sns.boxplot(y='Outcome', x='Jaccard', data=plotdf, orient="h", palette=sns.color_palette(cp1))
plt.xlabel('Jaccard Similarity')
plt.title('Morgan Fingerprint 8192 bits, Radius 3')
plt.tight_layout()
plt.show()

# clustermap
dfplot=pd.DataFrame(matsim, index=df['name'], columns=df['name'])
g = sns.clustermap(dfplot, row_colors=row_colors,col_colors=row_colors,
                   cmap="Greys",vmin=0, vmax=1,figsize=(50,50),
                   cbar_kws={"orientation": "horizontal", "label": "Jaccard"},
                   norm=colors.PowerNorm(gamma=0.5))
ax = g.ax_heatmap

# Draw the legend bar for the classes
for label in lut:
    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                            label=label, linewidth=0)
g.ax_col_dendrogram.legend( loc='right',ncol=1, fontsize=50,markerscale=2.)
# Adjust the postion of the main colorbar for the heatmap
g.cax.set_position([.35, 0.03, .5, .02])

'''dendro_box = g.ax_row_dendrogram.get_position()
dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
g.cax.set_position(dendro_box)'''
ax.set_xlabel('Compounds')
ax.set_ylabel('Compounds')
plt.tight_layout()
plt.savefig('./JaccardSimilarityexp_gen_corn.svg')

plt.show()


# plot a vclustermap of corn to kinetics keys
kincorndf=pd.read_csv('corn_top_kinetics.csv')

#need to get the smiles for all the keys. first get the keys
kcornkeys={}
kcorn=pd.DataFrame()
for item in kincorndf['inchi_keys'].values:
    item=item.replace('{','').replace('}','').replace('\\','').replace('"','')
    items=item.split(',')
    for ik in items:
        ik=ik.strip().replace("\'",'')
        ikl=ik.split('-')
        fb=ikl[0]
        kcornkeys[ik]=kincorndf[kincorndf['inchi_key']==fb].to_dict()

# adjust for reading a list from a text file
cornkeys=set()
with open('existing_fraction_corn.txt') as fin:
    for line in fin:
        line=line.strip()
        cornkeys.add(line)

for ik in mdm:
    if 'name' in mdm[ik]:
        if 'D-Fructose' in mdm[ik]['name']:
            print(ik)
            print(mdm[ik]['name'])
            if 'orgs' in mdm[ik] and 'zea mays' in mdm[ik]['orgs']:
                print('in corn')
        if 'L-mal' in mdm[ik]['name']:
            print(ik)
            print(mdm[ik]['name'])
            if 'orgs' in mdm[ik] and 'zea mays' in mdm[ik]['orgs']:
                print('in corn')
for kc in tqdm(cornkeys):
    kc=kc.strip()
    kc=kc.replace("\'",'')

    if kc in ndm['name'].keys():
        name=ndm['name'][kc]
    elif kc in mdm and 'name' in mdm[kc]:
        name=mdm[kc]['name']
    else:
        name=getName(kc)
    if kc in ik2isosmiles:
        smiles=ik2isosmiles[kc]
    else:
        smiles=getSmiles(kc)
    if kc in ik2type:
        if len(ik2type[kc])==2 or 'Secondary' in ik2type[kc]:
            cls = 'secondary'
        else:
            cls='primary'
    else:
        cls='unknown'
    '''for ke in kcornkeys[kc]['experiment'].keys():
        if kcornkeys[kc]['experiment'][ke]==1.0:
            exp=1
        else:
            exp=0'''
    exp=1
    values_to_add = {'name': name,'smiles':smiles,'inchi key':kc,'type':cls,'experiment':exp}
    row_to_add = pd.Series(values_to_add, name=kc)
    kcorn = kcorn.append(row_to_add)

kcorn.to_csv('kcorn_existing.csv')
kcorn.to_pickle('kcorn_existing.pkl')
# the input is ready - start analysis on the general list of keys
addMolecule(kcorn)
# remove all the lines that do not have a molecule
# split the dataframe into smaller chunks
getFpBitVec(kcorn)

kembeds=getEmbeddings(kcorn)
#mapk=createUMAP(kembeds)
matsimk,matdistk=getJacSimMtx(kembeds,kcorn)


classes=kcorn['type'].unique()
ordered_classes=sorted(classes)
counterST=collections.Counter(kcorn['type'])
counterSTkey=sorted(counterST, key=counterST.get, reverse=True)
counttypes=[counterST[k] for k in counterSTkey]
print([(k, counterST[k]) for k in counterSTkey])

lut = dict(zip(set(counterSTkey), sns.color_palette(cp1)))
row_colors = list(kcorn['type'].map(lut))

# clustermap
sns.set(font_scale=1)
dfplot1=pd.DataFrame(matsimk, index=kcorn['name'], columns=kcorn['name'])
g = sns.clustermap(dfplot1, row_colors=row_colors,col_colors=row_colors,
                   cmap="Greys",vmin=0, vmax=1,figsize=(50,50),
                   cbar_kws={"orientation": "horizontal", "label": "Jaccard"},
                   norm=colors.PowerNorm(gamma=0.5))
ax = g.ax_heatmap

# Draw the legend bar for the classes
for label in lut:
    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                            label=label, linewidth=0)
g.ax_col_dendrogram.legend( loc='right',ncol=1, fontsize=50,markerscale=2.)
# Adjust the postion of the main colorbar for the heatmap
g.cax.set_position([.35, 0.03, .5, .02])

'''dendro_box = g.ax_row_dendrogram.get_position()
dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
g.cax.set_position(dendro_box)'''
ax.set_xlabel('Compounds')
ax.set_ylabel('Compounds')
plt.tight_layout()
plt.savefig('./JaccardSimilarity_corn_existing.svg')

plt.show()