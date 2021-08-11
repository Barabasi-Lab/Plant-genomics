''' The new plant to try and evaluate the real significance of overlaps using PCA to count the effective compounds
in each subset and overall. The piepline is:
1. calculate the similarity matrix of the data subset
2. perform PCA
3. establish the number of components explaining the variance that is not random
4. get the effective number of compounds
5. calculate the overlap with the real original compound numbers
6. repeat items 1-4 for the overlap
7. plug the numbers into a hyhpergeometric test and get results'''


# imports
import sys
import scikitplot as skplt
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
from sklearn.decomposition import PCA
from sklearn import datasets
import math


# -------------------- functions -------------------
def convert2fb(iklist):
    ''' converts a full key list to first block list'''
    fblist=set()
    for ik in iklist:
        ikl=ik.split('-')
        fblist.add(ikl[0])
        fb2ik[ikl[0]]=ik
    return fblist

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
                if  inchikey in mdm and 'chebi' in mdm[inchikey]:
                    ccomp = ch.getLiteEntity(mdm[inchikey]['chebi'][0])
                    chcomp=ch.getCompleteEntity(mdm[inchikey]['chebi'][0])
                    if chcomp:
                        smiles=chcomp.smiles


                        mdm[inchikey]['smiles'] = smiles
                        break
                else:
                    if inchikey in ndm:
                        smiles=ndm[inchikey]['smiles']
                        break
                    else:
                        smiles=None
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

def getInputMtx(ikset):
    ''' adapted for first block'''
    ''' 1. prepare dataframe'''
    timestr=time.time()
    print(len(ikset))
    df=pd.DataFrame()
    for ik in tqdm(ikset):
        if ik in fb2ik:
            ik=fb2ik[ik]
        else:
            pass
        if ik in ik2isosmiles:
            smiles=ik2isosmiles[ik]
        else:
            smiles=getSmiles(ik)
        if  smiles and 'C1=2N3C(C=C4N5=C(C=C6N7C(=CC8=N(C(=C1)C(=C8CCC([O-])=O)C([H])=O)[Fe]735)C(=C6C)CCC([O-])=O)C(=C4C)C=C)=C(C2C)[C@H](CC/C=C(\C)/CC/C=C(/CCC=C(C)C)\C)O' in smiles:
            continue
        if ik in ik2name:
            name=ik2name[ik]
        else:
            name=ik
        if smiles:
            values_to_add = {'inchikey':ik, 'smiles': smiles,'name':name}
            row_to_add = pd.Series(values_to_add, name=org)
            df = df.append(row_to_add)
    addMolecule(df)
    getFpBitVec(df)
    genembeds = getEmbeddings(df)
    matsim, matdist = getJacSimMtx(genembeds, df)
    timeend=time.time()
    rt=timeend-timestr
    print(rt/60)
    return matsim


def calcPCA(mtx,title, ikset):
    ''' performs PCA, and save the significance plot'''
    pca = PCA()
    pca.fit_transform(mtx)
    '''skplt.decomposition.plot_pca_component_variance(pca)
    plt.show()'''
    # divide the eigenvalues by the sum of the number of compounds   (counter)
    deig=[]
    esum=0
    jsum=0
    effective_num=0 #pca.components_.shape[0]
    Ncomp=pca.components_.shape[0] #len(ikset)
    avg = 1 / len(ikset)
    sd = math.sqrt(((avg * (1 - avg)) / (len(ikset) + 1)))
    sd1=math.sqrt(((Ncomp-1)/(Ncomp+1))*(1/Ncomp**2))
    test = pca.explained_variance_ratio_
    for i,eig in enumerate(pca.explained_variance_ratio_):
        esum+=eig
        # james correction
        jsum+=(eig-1)**2
        deigt=eig/esum
        if deigt > avg:
            effective_num+=1
        deig.append(deigt)
    #jvar=jsum/(Ncomp-1)
    # james M effective
    #effective_num=round(Ncomp*(1-(Ncomp-1)*(jvar/Ncomp**2)),0)

    '''# plot
    var = np.cumsum(np.round(pca.explained_variance_, decimals=3) * 100)
    test=pca.explained_variance_ratio_
    test1=pca.components_
    plt.plot(deig[:40])
    plt.xlabel('number of components')
    plt.ylabel('Eigenvalues/sum(trace)')
    plt.title(title)
    plt.show()

    plt.plot(pca.explained_variance_ratio_)
    plt.xlabel('number of components')
    plt.ylabel('explained_variance_ratio')
    plt.title(title)
    plt.show()'''
    '''plt.rcParams["figure.figsize"] = (15, 6)

    fig, ax = plt.subplots()
    xi = np.arange(1, test.shape[0]+1, step=1)
    y = deig

    plt.ylim(0.0, 1.1)
    plt.plot(xi, y, marker='o', linestyle='--', color='b')

    plt.xlabel('Number of Components')
    plt.xticks(np.arange(0, test.shape[0], step=1))  # change from 0-based array index to 1-based human-readable label
    plt.ylabel('Eigenvalues/sum(trace)')
    plt.title('The number of components needed to explain variance')

    plt.axhline(y=0.95, color='r', linestyle='-')
    plt.text(0.5, 0.85, '95% cut-off threshold', color='red', fontsize=16)

    ax.grid(axis='x')
    plt.show()'''
    return deig[:15], effective_num




# ------------------- Main -------------------------
fb2ik={}
ch=ChEBI()
ndmdf=pd.read_csv('NDM_MasterJan21.csv')
ndm=ndmdf.set_index('InChIKey').to_dict()
mdm=pickle.load(open('mdm_040621.pkl','rb'))
ik2isosmiles=pickle.load(open('ik2isosmiles.pkl','rb'))
ik2name=pickle.load(open('ik2name_jun21.pkl','rb'))
orgs=pickle.load(open('orgkeys_new.pkl','rb'))
orgdgscores=pickle.load(open('kineticScores_full_key.pkl','rb'))
orgOpt=pickle.load(open('optimal_orgs.pkl','rb'))
orgdgscoresnew={}
for org in orgdgscores:
    orgl=org.lower()
    orgdgscoresnew[orgl]=orgdgscores[org]

orgdgscores=orgdgscoresnew

totonlyexpri=set()
plantlist=set()
# get big M for the overlaps with experiments
for org in orgs:
    if 'experiment' in orgs[org]:
        plantlist.add(org)
        for atr in orgs[org]:
            if atr in orgs[org]:
                totonlyexpri = totonlyexpri | orgs[org][atr]

totonlyexpri=convert2fb(totonlyexpri)
# step 1 - prepare input for similarity matrix calculations - specifically get smiles.
'''# first pass - estimate runtime for corn

corng=orgs['zea mays']['kegg']|orgs['zea mays']['pmn']
#manually curated inchi key discrepencies fix as a result of keys coming from different sources
cur = ['LFTYTUAZOPRMMI-MPIASZHXSA-N', 'SRBFZHDQGSBBOR-SOOFDHNKSA-N', 'MEFKEPWMEQBLKI-XCPQSEKJSA-O',
               'BAWFJGJZGIEFAR-NNYOXOHSSA-N', 'CBIDVWSRUUODHL-QTSLKERKSA-N', 'VOXXWSYKYCBWHO-UHFFFAOYSA-N',
               'JPIJQSOTBSSVTP-STHAYSLISA-M', 'OSJPPGNTCRNQQC-UHFFFAOYSA-N', 'NKDFYOWSKOHCCO-YPVLXUMRSA-N',
               'JZNWSCPGTDBMEW-UHFFFAOYSA-N', 'SBLKVIQSIHEQOF-UPHRSURJSA-N', 'RBNPOMFGQQGHHO-UHFFFAOYSA-N',
               'JYJIGFIDKWBXDU-MNNPPOADSA-N', 'FYSSBMZUBSBFJL-UHFFFAOYSA-N', 'ACWASDPGAVYCNI-JKUQZMGJSA-N',
               'BRGMHAYQAZFZDJ-PVFLNQBWSA-N', 'LLCSXHMJULHSJN-UHFFFAOYSA-N', 'ODBLHEXUDAPZAU-ZAFYKAAXSA-N',
               'HEBKCHPVOIAQTA-NGQZWQHPSA-N', 'RHGKLRLOHDJJDR-UHFFFAOYSA-N', 'HXEACLLIILLPRG-UHFFFAOYSA-N',
               'ZZLHPCSGGOGHFW-ZMQIUWNVSA-N', 'VCWMRQDBPZKXKG-SPBUTQSFSA-N', 'JFCQEDHGNNZCLN-UHFFFAOYSA-N',
               'BJEPYKJPYRNKOW-UHFFFAOYSA-N', 'LKDRXBCSQODPBY-VRPWFDPXSA-N']
corng=corng|set(cur)

corndbs=orgs['zea mays']['foodb']|orgs['zea mays']['dfc']|orgs['zea mays']['usda']
cornexp=orgs['zea mays']['experiment']
cornkinetics=set()
for dgik in orgdgscores['zea mays']:
    if orgdgscores['zea mays'][dgik]['sum']>0:
        cornkinetics.add(dgik)
corntotal=corng | corndbs | cornexp | cornkinetics

matcorng=getInputMtx(corng)
lg,gnum=calcPCA(matcorng,'Genome')
matcdbs=getInputMtx(corndbs)
ldbs, dbsnum=calcPCA(matcdbs,'Databases')
matexp=getInputMtx(cornexp)
lexp, expnum=calcPCA(matexp,'Experiment')
matkinetics=getInputMtx(cornkinetics)
lk,knum=calcPCA(matkinetics,'Kinetics')
matcorntot=getInputMtx(corntotal)
ltot,totnum=calcPCA(matcorntot,'Total')

# plot all the pca outputs
p1=plt.plot(lg, color='green', label='Genome')
p2=plt.plot(ldbs, color='blue', label='Databases')
p3=plt.plot(lexp, color='orange', label='Experiments')
p4=plt.plot(lk, color='grey', label='kinetics')
p5=plt.plot(ltot, color='black', label='total')
plt.xlabel('number of components')
plt.ylabel('Eigenvalues/sum(trace)')
plt.legend()

plt.show()'''

# big M
matbigm=getInputMtx(totonlyexpri)
lbigm,nbigm=calcPCA(matbigm,'Big M',totonlyexpri)

#nbigm=168



# expand to all the plants and plot
resdf=pd.DataFrame()
for org in plantlist:
    if 'pyrus communis' in org:
        continue
    # genome set
    if 'kegg' not in orgs[org]:
        gset=orgs[org]['pmn']
    if 'pmn' not in orgs[org]:
        gset = orgs[org]['kegg']
    if 'kegg' in orgs[org] and 'pmn' in orgs[org]:
        gset = orgs[org]['kegg'] | orgs[org]['pmn']
    # manually curated inchi key discrepencies with NDM fix as a result of keys describing the same molecule but coming from different sources (unichem, PMN)
    cur = ['LFTYTUAZOPRMMI-MPIASZHXSA-N', 'SRBFZHDQGSBBOR-SOOFDHNKSA-N', 'MEFKEPWMEQBLKI-XCPQSEKJSA-O',
           'BAWFJGJZGIEFAR-NNYOXOHSSA-N', 'CBIDVWSRUUODHL-QTSLKERKSA-N', 'VOXXWSYKYCBWHO-UHFFFAOYSA-N',
           'JPIJQSOTBSSVTP-STHAYSLISA-M', 'OSJPPGNTCRNQQC-UHFFFAOYSA-N', 'NKDFYOWSKOHCCO-YPVLXUMRSA-N',
           'JZNWSCPGTDBMEW-UHFFFAOYSA-N', 'SBLKVIQSIHEQOF-UPHRSURJSA-N', 'RBNPOMFGQQGHHO-UHFFFAOYSA-N',
           'JYJIGFIDKWBXDU-MNNPPOADSA-N', 'FYSSBMZUBSBFJL-UHFFFAOYSA-N', 'ACWASDPGAVYCNI-JKUQZMGJSA-N',
           'BRGMHAYQAZFZDJ-PVFLNQBWSA-N', 'LLCSXHMJULHSJN-UHFFFAOYSA-N', 'ODBLHEXUDAPZAU-ZAFYKAAXSA-N',
           'HEBKCHPVOIAQTA-NGQZWQHPSA-N', 'RHGKLRLOHDJJDR-UHFFFAOYSA-N', 'HXEACLLIILLPRG-UHFFFAOYSA-N',
           'ZZLHPCSGGOGHFW-ZMQIUWNVSA-N', 'VCWMRQDBPZKXKG-SPBUTQSFSA-N', 'JFCQEDHGNNZCLN-UHFFFAOYSA-N',
           'BJEPYKJPYRNKOW-UHFFFAOYSA-N', 'LKDRXBCSQODPBY-VRPWFDPXSA-N','HSCJRCZFDFQWRP-ABVWGUQPSA-N',
           'AAWZDTNXLSGCEK-LNVDRNJUSA-N','ILGMGHZPXRDCCS-UHFFFAOYSA-N','NGSWKAQJJWESNS-UHFFFAOYSA-N',
           'FEWJPZIEWOKRBE-UHFFFAOYSA-N','PADQINQHPQKXNL-UHFFFAOYSA-N','FDHFJXKRMIVNCQ-OSIZZBRKSA-N',
           'QQHJDPROMQRDLA-UHFFFAOYSA-N','QAIPRVGONGVQAS-RQOWECAXSA-N','YDBYJHTYSHBBAU-UHFFFAOYSA-O',
           'POJWUDADGALRAB-UHFFFAOYSA-N','CUZKLRTTYZOCSD-UHFFFAOYSA-N','BJRNKVDFDLYUGJ-UHFFFAOYSA-N']
    gset = gset | set(cur)
    # convert to first block
    gset=convert2fb(gset)


    # databases set
    dbset=orgs[org]['foodb'] | orgs[org]['dfc'] | orgs[org]['usda']
    # convert to first block
    dbset=convert2fb(dbset)


    # experiment set
    expset=orgs[org]['experiment']
    expset=convert2fb(expset)
    # kinetics set
    kset=set()
    for dgik in orgdgscores[org]:
        if orgdgscores[org][dgik]['average'] > 0:
            kset.add(dgik)
    kset=convert2fb(kset)
    # kinetics optimal fraction
    okset=orgOpt[org]
    # convert to fb
    okset=convert2fb(okset)

    totalset = gset | dbset | expset | kset
    matg = getInputMtx(gset)
    lg, gnum = calcPCA(matg, 'Genome',gset)
    matdbs = getInputMtx(dbset)
    ldbs, dbsnum = calcPCA(matdbs, 'Databases',dbset)
    matexp = getInputMtx(expset)
    lexp, expnum = calcPCA(matexp, 'Experiment',expset)
    matkinetics = getInputMtx(kset)
    lk, knum = calcPCA(matkinetics, 'Kinetics',kset)
    mattot = getInputMtx(totalset)
    ltot, totnum = calcPCA(mattot, 'Total',totalset)

    # overlaps
    # genome vs. experiments
    gve = gset & expset
    matgve = getInputMtx(gve)
    gvematu=np.triu(matgve)
    gvematu.flatten()
    gvematlist=np.ndarray.tolist(gvematu)
    gvematlist = [i for i in gvematlist if i != 0]
    gvematlist = [i for i in gvematlist if i != 1.0]
    # plot distribution
    sns.distplot(gvematlist)
    plt.show()
    lgve, gvenum = calcPCA(matgve, 'Genome vs. exp',gve)


    # test existing
    pvalge = hypergeom.sf(gvenum - 1, nbigm, expnum, gnum)
    if pvalge>0:
        pvalge=abs(math.log10(pvalge))
    expectedge = (expnum * gnum) / nbigm
    prbge = hypergeom.cdf(gvenum, nbigm,expnum, gnum)
    gveratio=gvenum/expectedge


    # real numbers genome vs. experiments
    rpvalge = hypergeom.sf(len(gve) - 1, len(totonlyexpri), len(expset), len(gset))
    rexpectedge = (len(expset) * len(gset)) / len(totonlyexpri)
    rprbge = hypergeom.cdf(len(gve), len(totonlyexpri), len(expset), len(gset))
    rgveratio = len(gve) / rexpectedge


    # databases vs. experiments
    dbve = dbset & expset
    matdbve = getInputMtx(dbve)
    ldbve, dbvenum = calcPCA(matdbve, 'databases vs.exp',dbve)

    # test existing
    pvaldbe = hypergeom.sf(dbvenum - 1, nbigm, expnum, dbsnum)
    if pvaldbe>0:
        pvaldbe=abs(math.log10(pvaldbe))
    expecteddbe = (expnum * dbsnum) / nbigm
    prbdbe = hypergeom.cdf(dbvenum, nbigm, expnum, dbsnum)
    dbveratio=dbvenum/expecteddbe


    # real numbers dbs vs. experiments
    rpvaldbe = hypergeom.sf(len(dbve) - 1, len(totonlyexpri), len(expset), len(dbset))
    rexpecteddbe = (len(expset) * len(dbset)) / len(totonlyexpri)
    rprbdbe = hypergeom.cdf(len(dbve), len(totonlyexpri), len(expset), len(dbset))
    rdbveratio = len(dbve) / rexpecteddbe


    # kinetics vs. experiments
    kve = kset & expset
    matkve = getInputMtx(kve)
    lkve, kvenum = calcPCA(matkve, 'kinetics vs. exp',kve)

    # test existing
    pvalke = hypergeom.sf(kvenum - 1, nbigm, expnum, knum)
    if pvalke>0:
        pvalke=abs(math.log10(pvalke))
    expectedke = (expnum * knum) / nbigm
    prbke = hypergeom.cdf(kvenum, nbigm, expnum, knum)
    kveratio=kvenum/expectedke

    # real numbers kinetics vs. experiments
    rpvalke = hypergeom.sf(len(kve) - 1, len(totonlyexpri), len(expset), len(kset))
    rexpectedke = (len(expset) * len(kset)) / len(totonlyexpri)
    rprbke = hypergeom.cdf(len(kve), len(totonlyexpri), len(expset), len(kset))
    rkveratio = len(kve) / rexpectedke



    # databases vs. genome
    dbvg = gset & dbset
    matdbvg = getInputMtx(dbvg)
    ldbvg, dbvgnum = calcPCA(matdbvg, 'databases vs. genome',dbvg)

    # test existing
    pvaldbvg = hypergeom.sf(dbvgnum - 1, nbigm, dbsnum, gnum)
    expecteddbvg = (dbsnum * gnum) / nbigm
    prbdbvg = hypergeom.cdf(dbvgnum, nbigm, dbsnum, gnum)
    dbvgratio=dbvgnum/expecteddbvg

    # real numbers dbs vs. genome
    rpvaldbvg = hypergeom.sf(len(dbvg) - 1, len(totonlyexpri), len(dbset), len(gset))
    rexpecteddbvg = (len(dbset) * len(gset)) / len(totonlyexpri)
    rprbdbvg = hypergeom.cdf(len(dbvg), len(totonlyexpri), len(dbset), len(gset))
    rdbvgratio = len(dbvg) / expecteddbvg


    # databases vs. kinetics
    dbvk = kset & dbset
    matdbvk = getInputMtx(dbvk)
    ldbvk, dbvknum = calcPCA(matdbvk, 'kinetics vs. databases',dbvk)

    # test existing
    pvaldbvk = hypergeom.sf(dbvknum - 1, nbigm, dbsnum, knum)
    expectedbvk = (dbsnum * knum) / nbigm
    prbdbvk = hypergeom.cdf(dbvknum, nbigm, dbsnum, knum)
    dbvkratio=dbvknum/expectedbvk

    # real numbers dbs vs. kinetics
    rpvaldbvk = hypergeom.sf(len(dbvk) - 1, len(totonlyexpri), len(dbset), len(kset))
    rexpectedbvk = (len(dbset) * len(kset)) / len(totonlyexpri)
    rprbdbvk = hypergeom.cdf(len(dbvk), len(totonlyexpri), len(dbset), len(kset))
    rdbvkratio = len(dbvk) / rexpectedbvk




    # update the dataframe
    values_to_add = {'plant': org, 'M':nbigm,'effective number genome':gnum,'effective number dbs':dbsnum,'effective number kinetics':knum,
                     'effective number experiments': expnum,
                     'Real overlap genome vs. experiments': len(gve), 'Effective overlap genome vs. experiments': gvenum,
                     'P-value genome vs. experiments': pvalge,'cdf genome vs. experiment':prbge,
                     'Real overlap databases vs. experiments': len(dbve),
                     'Effective overlap databases vs. experiments': dbvenum, 'P-value databases vs. experiments': pvaldbe,
                     'cdf databases vs. experiments': prbdbe,
                     'Real overlap kinetics vs. experiments': len(kve), 'Effective overlap kinetics vs. experiment': kvenum,
                     'P-value kinetics vs. experiments': pvalke, 'cdf kinetics vs. experiments':prbke,
                     'Real overlap databases vs. genome': len(dbvg),
                     'Effective overlap databases vs. genome': dbvgnum,
                     'P-value databases vs. genome': pvaldbvg,'cdf databases vs. genome': prbdbvg,
                     'Real overlap kinetics vs. databases': len(dbvk),
                     'Effective overlap kinetics vs. genome': dbvknum,
                     'P-value kinetics vs. databases': pvaldbvk, 'cdf kinetics vs. databases':prbdbvk,
                     'genome vs. experiments ratio': gveratio,'databases vs. experiments ratio': dbveratio,
                     'kinetics vs. experiments ratio': kveratio, 'genome vs. databases ratio':dbvgratio,
                     'kinetics vs. databases ratio': dbvkratio,
                     'genome fraction':len(gset), 'databases fraction':len(dbset),'kinetics fraction':len(kset),
                     'experiments fraction': len(expset),'real big M': len(totonlyexpri),
                     'real P-value genomes vs. experiments': rpvalge,'real genome vs. experiments ratio':rgveratio,
                     'real P-value databases vs. experiments': rpvaldbe, 'real databases vs. experiments ratio': rdbveratio,
                     'real P-value kinetics vs. experiments': rpvalke, 'real kinetics vs. experiments ratio': rkveratio,
                     'real P-value databases vs. genome': rpvaldbvg,'real genome vs. databases ratio': rdbvgratio,
                     'real P-value kinetics vs. databases': rpvaldbvk, 'real kinetics vs. databases ratio': rdbvkratio
                     }
    row_to_add = pd.Series(values_to_add, name=org)
    resdf = resdf.append(row_to_add)


# plot p-values
sns.set(font_scale=3)
plt.figure(figsize=(20, 20))
f1=sns.scatterplot(data=resdf, x="plant", y="P-value genome vs. experiments", color='steelblue', s=250, label='genome vs. experiments')
f2=sns.scatterplot(data=resdf, x="plant", y="P-value databases vs. experiments", color='yellowgreen', s=250, label='databases vs. experiments')
f3=sns.scatterplot(data=resdf, x="plant", y="P-value kinetics vs. experiments", color='goldenrod',s=250,label='kinetics vs. experiments')
'''f4=sns.scatterplot(data=resdf, x="plant", y="P-value kinetics vs. databases", color='blueviolet',s=200, label='kinetics vs. databases')
f5=sns.scatterplot(data=resdf, x="plant", y="P-value databases vs. genome", color='red',s=200,label='databases vs. genome')'''
'''f1.set(yscale="log")
f2.set(yscale="log")
f3.set(yscale="log")'''
plt.axhline(y=abs(math.log10(0.05)), color='black', linestyle='--',linewidth=5)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=3, mode="expand", borderaxespad=0., fontsize=26)
plt.draw()
f1.set_xticklabels(f1.get_xticklabels(), rotation=90)

plt.tight_layout()
plt.savefig('pca_p_vals_james.svg')
plt.show()


# plot ratios
sns.set(font_scale=3)
plt.figure(figsize=(20, 20))
f1=sns.scatterplot(data=resdf, x="plant", y="genome vs. experiments ratio", color='steelblue', s=250, label='genome vs. experiments')
f2=sns.scatterplot(data=resdf, x="plant", y="databases vs. experiments ratio", color='yellowgreen', s=250, label='databases vs. experiments')
f3=sns.scatterplot(data=resdf, x="plant", y="kinetics vs. experiments ratio", color='goldenrod',s=250,label='kinetics vs. experiments"')
'''f4=sns.scatterplot(data=resdf, x="plant", y="kinetics vs. databases ratio", color='blueviolet',s=200, label='kinetics vs. databases')
f5=sns.scatterplot(data=resdf, x="plant", y="genome vs. databases ratio", color='red',s=200,label='databases vs. genome')'''
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=3, mode="expand", borderaxespad=0., fontsize=26)
plt.draw()
f1.set_xticklabels(f1.get_xticklabels(), rotation=90)
plt.axhline(y=1, color='black', linestyle='--',linewidth=5)
plt.tight_layout()
plt.savefig('pca_ratio_james.svg')
plt.show()


# plot p-values
sns.set(font_scale=3)
plt.figure(figsize=(20, 20))
f1=sns.scatterplot(data=resdf, x="plant", y="real P-value genomes vs. experiments", color='steelblue', s=250, label='genome vs. experiments')
f2=sns.scatterplot(data=resdf, x="plant", y="real P-value databases vs. experiments", color='yellowgreen', s=250, label='databases vs. experiments')
f3=sns.scatterplot(data=resdf, x="plant", y="real P-value kinetics vs. experiments", color='goldenrod',s=250,label='kinetics vs. experiments')
'''f4=sns.scatterplot(data=resdf, x="plant", y="real P-value kinetics vs. databases", color='blueviolet',s=200, label='kinetics vs. databases')
f5=sns.scatterplot(data=resdf, x="plant", y="real P-value databases vs. genome", color='red',s=200,label='databases vs. genome')'''
f1.set(yscale="log")
plt.axhline(y=0.05, color='black', linestyle='--',linewidth=5)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=3, mode="expand", borderaxespad=0., fontsize=26)
plt.draw()
f1.set_xticklabels(f1.get_xticklabels(), rotation=90)

plt.tight_layout()
plt.savefig('pca_p_vals_real_f.svg')
plt.show()


# plot ratios
sns.set(font_scale=3)
plt.figure(figsize=(25, 25))
f1=sns.scatterplot(data=resdf, x="plant", y="real genome vs. experiments ratio", color='steelblue', s=250, label='genome vs. experiments')
f2=sns.scatterplot(data=resdf, x="plant", y="real databases vs. experiments ratio", color='yellowgreen', s=250, label='databases vs. experiments')
f3=sns.scatterplot(data=resdf, x="plant", y="real kinetics vs. experiments ratio", color='goldenrod',s=250,label='kinetics vs. experiments"')
'''f4=sns.scatterplot(data=resdf, x="plant", y="real kinetics vs. databases ratio", color='blueviolet',s=200, label='kinetics vs. databases')
f5=sns.scatterplot(data=resdf, x="plant", y="real genome vs. databases ratio", color='red',s=200,label='databases vs. genome')'''
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=3, mode="expand", borderaxespad=0., fontsize=26)
plt.draw()
f1.set_xticklabels(f1.get_xticklabels(), rotation=90)
plt.axhline(y=1, color='black', linestyle='--',linewidth=5)
plt.tight_layout()
plt.savefig('pca_ratio_real_f.svg')
plt.show()



print('done')


