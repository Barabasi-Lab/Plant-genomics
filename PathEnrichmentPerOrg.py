'''Profiling the coverage of KEGG pathways- how many compounds are in each pathway.
Doing a hypergeometric test to test for enrichment (biases for certain pathways.
That would be the baseline. Step 2 - collecting KEGG annotations from MDM and seeing that I have more than what
was found in KEGG.
Recounting pathway coverages - how many compounds per pathway.
Redoing hypergeometric test to test for a shift

Written by Shany Ofaim, the Barabasi lab CCNR Northeastern University 2021'''


# imports
import plotly.express as px
import math
import pickle
import re
import numpy as np
from varname import nameof
from bioservices import *
import json
from scipy.stats import hypergeom
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.sandbox.stats.multicomp import multipletests
# ----- functions ---------
def cleanhtml(raw_html):
    # html tag cleaner copied from stack overflow, cleans non tag entities as well
  cleanr = re.compile('<.*?>|&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-f]{1,6});')
  cleantext = re.sub(cleanr, '', raw_html)
  return cleantext
def plotlyPCA(df):
    fig = px.scatter(df, x=df['PC1'], y=df['PC2'], color=df['cat'])
    fig.show()

def plotlyPCA3d(df):
    fig = px.scatter_3d(
        df, x=df['PC1'], y=df['PC2'], z=df['PC3'], color=df['plant']
    )
    fig.update_traces(marker=dict(size=6,
                                  line=dict(width=1,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
    fig.show()


def plotPCA(dframe,fname):

    # plots a pca plot for a dataframe using a scaled version of the data
    plants = dframe.plant.tolist()
    paths = dframe.pathway.tolist()
    cats = dframe.cat.tolist()
    types = dframe.cls.tolist()
    cats_set=set(cats)
    plants_set=set(plants)
    orgs_data = dframe.select_dtypes(np.number)
    # scale the data
    random_state = 0
    pca_scaled = make_pipeline(StandardScaler(),
                               PCA(n_components=4, random_state=random_state))
    orgs_pca_scaled = pca_scaled.fit_transform(orgs_data)
    pc_df_scaled = pd.DataFrame(data=orgs_pca_scaled,
                                columns=['PC1', 'PC2', 'PC3', 'PC4'])
    pc_df_scaled['plant'] = plants
    pc_df_scaled['pathway'] = paths
    pc_df_scaled['cat'] = cats
    pc_df_scaled['types'] = types
    paths_set = set(paths)
    explained_scaled = pca_scaled.named_steps['pca'].explained_variance_ratio_ * 100

    plt.figure(figsize=(30, 20))
    sns.set(font_scale=2.5)
    sns.scatterplot(x="PC1", y="PC2",
                    data=pc_df_scaled,
                    hue="cat", palette=sns.color_palette('viridis', n_colors=len(cats_set)), s=300, style='types',
                    style_order=['Primary','Secondary'])
    plt.xlabel("PC1: " + f'{explained_scaled[0]:.0f}' + "%")
    plt.ylabel("PC2: " + f'{explained_scaled[1]:.0f}' + "%")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=5, mode="expand", borderaxespad=0., fontsize=16)
    plt.tight_layout()
    plt.show()
    filname=fname+'.csv'
    pc_df_scaled.to_csv(filname)
    plotlyPCA(pc_df_scaled)
    #plotlyPCA3d(pc_df_scaled)

def getPathways():
    pass
def cleanhtml(raw_html):
    # html tag cleaner copied from stack overflow, cleans non tag entities as well
  cleanr = re.compile('<.*?>|&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-f]{1,6});')
  cleantext = re.sub(cleanr, '', raw_html)
  return cleantext
# ----------- main -----------
kegg_con=KEGG()
# read the pathway and reaction mapping for all organisms
cat2class={}
# primary secondary mapping
with open('Primary_Secondary.txt') as fin:
    for line in fin:
        line=line.split('\t')
        cat=line[0].strip()
        clas=line[1].strip()
        cat2class[cat]=clas
# read mdm
print('loading MDM')
mdm={}
with open('MDM_021721.json') as mdmin:
    mdm=json.load(mdmin)

kpath2code={}
specialCase={}
with open('pathwayCodes.txt') as fin:
    for line in fin:
        line=line.split('\t')
        code=line[0].strip().zfill(5)
        path=line[1].strip()
        if ',' in path:
            specialCase[path]=path.replace(',','')
        kpath2code[path]=code

# get the general pathway population for PMN
pmn_M=set()
genPathPMN={}
with open('Copy_of_All_pathways_of_PlantCyc.txt') as pmnin:
    for line in tqdm(pmnin):
        if line.startswith('Pathways'):
            continue
        line=line.split('\t')
        line[0]=line[0].replace('&','').replace(';','')
        path=cleanhtml(line[0])
        path=path.strip()
        pcode=line[5].strip()
        compounds=line[6].strip().split('//')
        compset=set()
        for comp in compounds:
            comp=comp.strip()
            if '' is  comp:
                continue
            comp=cleanhtml(comp)
            compset.add(comp.lower())
        pmn_M=pmn_M|(compset)
        if path not in genPathPMN:
            path=path.replace('"','')
            path=cleanhtml(path)
            genPathPMN[path]={}
        if 'compounds' not in genPathPMN[path]:
            genPathPMN[path]['compounds']=set()
        genPathPMN[path]['compounds']=genPathPMN[path]['compounds']|compset
        genPathPMN[path]['id']=pcode


# create the complete list of kegg pathways and a list of pmn pathways from mdm
pathway2cat={}
keggpathways=set()
pmnpathways=set()
org2pathways={}
mdmpath2comp={}
cnums=set()
cnums_total=set()
pnums=set()
pnums_total=set()
path2class={}
ik2pathways={}
ik2type={}
for ik in tqdm(mdm):
    if 'kegg' in mdm[ik]:
        cnums_total.add(ik)
    if 'plantcyc' in mdm[ik]:
        pnums_total.add(ik)
    if 'orgs' in mdm[ik]:
        for org in mdm[ik]['orgs']:
            # count c numbers
            if 'kegg' in mdm[ik]:
                cnums.add(ik)
            if 'plantcyc' in mdm[ik]:
                pnums.add(ik)
            if org not in org2pathways:
                org2pathways[org]={}
            for r in mdm[ik]['orgs'][org]:
                # get pathways
                pa = mdm[ik]['orgs'][org][r]
                for p in pa:
                    if not isinstance(pa[p], dict):
                        # pmn
                        if 'no mapping' in p:
                            continue
                        if 'pmn' not in org2pathways[org]:
                            org2pathways[org]['pmn']={}
                        if 'pmn' not in mdmpath2comp:
                            mdmpath2comp['pmn']={}
                        if p not in mdmpath2comp['pmn']:
                            mdmpath2comp['pmn'][p]=set()
                        if p not in org2pathways[org]['pmn']:
                            org2pathways[org]['pmn'][p]={}
                        if 'compounds' not in  org2pathways[org]['pmn'][p]:
                            org2pathways[org]['pmn'][p]['compounds']=set()
                        if 'iks' not in org2pathways[org]['pmn'][p]:
                            org2pathways[org]['pmn'][p]['iks']=set()
                        pmncomp=mdm[ik]['plantcyc'].lower()
                        org2pathways[org]['pmn'][p]['compounds'].add(pmncomp)
                        org2pathways[org]['pmn'][p]['iks'].add(ik)
                        mdmpath2comp['pmn'][p].add(mdm[ik]['plantcyc'].lower())
                        path2class[p]=cat2class[pa[p]]
                        if ik not in ik2type:
                            ik2type[ik]=set()
                        ik2type[ik].add(cat2class[pa[p]])

                        if p in genPathPMN and len(mdmpath2comp['pmn'][p])> len(genPathPMN[p]['compounds']):
                            pass
                        pmnpathways.add(p)
                        pathway2cat[p]=pa[p]
                        if ik not in ik2pathways:
                            ik2pathways[ik]=set()
                        ik2pathways[ik].add(p)
                    else:
                        # a kegg situation, r is infact an enzyme
                        for path in pa[p]:
                            if 'no mapping' in path:
                                continue
                            if 'kegg' not in org2pathways[org]:
                                org2pathways[org]['kegg']={}
                            if 'kegg' not in mdmpath2comp:
                                mdmpath2comp['kegg']={}
                            if path not in mdmpath2comp['kegg']:
                                mdmpath2comp['kegg'][path]=set()
                            if path not in org2pathways[org]['kegg']:
                                org2pathways[org]['kegg'][path]={}
                            if 'compounds' not in org2pathways[org]['kegg'][path]:
                                org2pathways[org]['kegg'][path]['compounds']=set()
                            if 'iks' not in org2pathways[org]['kegg'][path]:
                                org2pathways[org]['kegg'][path]['iks']=set()
                            if isinstance(mdm[ik]['kegg'],list):
                                mdm[ik]['kegg']=mdm[ik]['kegg'][0]
                            org2pathways[org]['kegg'][path]['compounds'].add(mdm[ik]['kegg'])
                            org2pathways[org]['kegg'][path]['iks'].add(ik)
                            mdmpath2comp['kegg'][path].add(mdm[ik]['kegg'])
                            keggpathways.add(path)
                            if ik not in ik2pathways:
                                ik2pathways[ik] = set()
                            ik2pathways[ik].add(path)
                            pathway2cat[path]=pa[p][path]
                            path2class[path]=cat2class[pa[p][path]]
                            if ik not in ik2type:
                                ik2type[ik] = set()
                            ik2type[ik].add(cat2class[pa[p][path]])

pickle.dump( ik2type, open( "ik2prim_sec.pkl", "wb" ) )
pickle.dump( ik2pathways, open( "ik2pathways.pkl", "wb" ) )
pickle.dump( path2class, open( "path2class.pkl", "wb" ) )
# get the general pathway population for KEGG
kegg_M=set()
genPathKEGG=pickle.load(open('genpaths_kegg.pkl','rb'))

for path in genPathKEGG:
    kegg_M=kegg_M|genPathKEGG[path]['compounds']
'''for p in tqdm(keggpathways):
    p=p.replace('"','')
    pa='rn'+kpath2code[p]
    path = kegg_con.parse(kegg_con.get(pa))
    # break down the reactions to compounds so the cofactors will be included.
    if 'REACTION' in path:
        compounds=set()
        for reac in path['REACTION']:
            if not reac.startswith('R'):
                continue
            kreac=kegg_con.parse(kegg_con.get(reac))
            # collect all the compounds for the reaction
            if 'EQUATION' in kreac:
                equation=kreac['EQUATION'].split(' <=> ')
                for eq in equation:
                    eq=eq.split(' + ')
                    for ecomp in eq:
                        ecomp=ecomp.strip()
                        compounds.add(ecomp)
        if compounds:
            if p not in genPathKEGG:
                genPathKEGG[p] = {}
            kegg_M = kegg_M | compounds
            genPathKEGG[p]['compounds'] = compounds
            genPathKEGG[p]['name'] = path['NAME']
            cat = path['CLASS'].split(';')
            genPathKEGG[p]['class'] = cat[1].strip()'''







# create the data for compound coverage in pathways for each organisms and calculate the hypergeopmetric test results
# collect the results in a dataframe and visualize

#pickle.dump( genPathKEGG, open( "genpaths_kegg.pkl", "wb" ) )

print('Total c numbers in mdm: '+ str(len(cnums_total)))
print('Total c numbers with orgs in mdm: '+ str(len(cnums)))
genpaths={}
genpaths.update(genPathPMN)
genpaths.update(genPathKEGG)
pickle.dump( genpaths, open( "genpaths.pkl", "wb" ) )
genEnriched={}
# Second pass - calculate hypergeometric p-value for each pathway - all organisms
genEnriched['kegg']={}
genEnriched['pmn']={}
genpath2fold={}
#
for pt in genPathKEGG:

    pts=pt
    N=n=len(genPathKEGG[pts]['compounds'])
    if pt in specialCase:
        pts=pt
        pt=specialCase[pt]
    # check the pathway is in the database
    if pt not in mdmpath2comp['kegg'] or pts not in mdmpath2comp['kegg']:
        continue
    #n=len(mdmpath2comp['kegg'][pt])
    if n>N:
        continue
    x=len(mdmpath2comp['kegg'][pt] & genPathKEGG[pts]['compounds'])
    pval = hypergeom.sf(x - 1, len(kegg_M), n, N)
    expected=(n*N)/len(kegg_M)
    fold=x/N
    genEnriched['kegg'][pt]={}
    genEnriched['kegg'][pt]=np.log10(pval)
    genpath2fold[pt]=fold
for pt in genPathPMN:
    # check the pathway exists in the database
    if pt not in mdmpath2comp['pmn']:
        continue
    N=n=len(genPathPMN[pt]['compounds'])
    #n=len(mdmpath2comp['pmn'][pt])
    if n>N:
        # plot the venn to investigate
        diff=mdmpath2comp['pmn'][pt]-genPathPMN[pt]['compounds']


    x=len(mdmpath2comp['pmn'][pt] & genPathPMN[pt]['compounds'])
    pval = hypergeom.sf(x - 1, len(pmn_M), n, N)
    expected = (n * N) / len(kegg_M)
    fold = x / N
    genEnriched['pmn'][pt]={}
    genEnriched['pmn'][pt]=np.log10(pval)
    genpath2fold[pt]=fold
# get the distribution of the enrichment for the general part
# turn into a dataframe
gendfkegg=pd.DataFrame(genEnriched['kegg'].items(), columns=['Pathway','Log10(P-value)'])
gendfpmn=pd.DataFrame(genEnriched['pmn'].items(), columns=['Pathway','Log10(P-value)'])
# plot the distribution
#distribution
plt.figure(figsize=(12, 10))
sns.set(font_scale=2.5)
g=sns.distplot(gendfkegg['Log10(P-value)'], color='dodgerblue',kde_kws=dict(linewidth=5), label='KEGG')
g1=sns.distplot(gendfpmn['Log10(P-value)'], color='red',kde_kws=dict(linewidth=5), label='PlantCyc')
plt.tight_layout()
plt.legend()
plt.show()




# third pass - getting the enriched pathways for each organism
org2enriched={}
org2numenriched={}
org2pathfold={}
totpaths=set()
ep2org={}
bad=set()
for org in tqdm(org2pathways):
    # kegg part
    if 'kegg' in org2pathways[org]:

        for pt in org2pathways[org]['kegg']:
            if any(x in pt for x in['Prodigiosin biosynthesis','Biosynthesis of vancomycin group antibiotics',
                                    'Acarbose and validamycin biosynthesis','engineered','Monobactam','monobactam']):
                continue
            totpaths.add(pt)
            pts = pt
            if pts in genPathKEGG:
                N =n= len(genPathKEGG[pts]['compounds'])
                if pt in specialCase:
                    pts = pt
                    pt = specialCase[pt]
                # check the pathway is in the database
                if pt not in org2pathways[org]['kegg'] or pts not in org2pathways[org]['kegg']:
                    continue
                #n = len(org2pathways[org]['kegg'][pt])
                if n > N:
                    continue
                x = len(set(org2pathways[org]['kegg'][pt]['compounds']) & genPathKEGG[pts]['compounds'])
                pval = hypergeom.sf(x - 1, len(kegg_M), n, N)
                expected = (n * N) / len(kegg_M)
                fold = x / N
                if org not in org2enriched:
                    org2enriched[org]={}
                org2enriched[org][pt] = {}
                org2enriched[org][pt]['pval']=pval
                org2enriched[org][pt]['coverage']=fold
                org2enriched[org][pt]['total in pathway']=N
                org2enriched[org][pt]['no. sampled']=x
                org2enriched[org][pt]['expected']=expected
                org2enriched[org][pt]['set']='KEGG'
                org2enriched[org][pt]['cat']=pathway2cat[pt]
                org2enriched[org][pt]['cls']=path2class[pt]
                if pval<0.5:
                    if org not in org2numenriched:
                        org2numenriched[org]=set()
                    org2numenriched[org].add(pt)
                    if pt not in ep2org:
                        ep2org[pt]=set()
                    ep2org[pt].add(org)

    # pmn part
    if 'pmn' in org2pathways[org]:

        for pt in org2pathways[org]['pmn']:
            if any(x in pt for x in ['engineered']):
                continue
            totpaths.add(pt)
            # check the pathway exists in the database
            if pt not in org2pathways[org]['pmn']:
                continue
            if pt in genPathPMN:

                N =n= len(genPathPMN[pt]['compounds'])
                #n = len(org2pathways[org]['pmn'][pt])
                if n > N:
                    # plot the venn to investigate
                    diff = org2pathways[org]['pmn'][pt]['compounds'] - genPathPMN[pt]['compounds']

                x = len(set(org2pathways[org]['pmn'][pt]['compounds']) & genPathPMN[pt]['compounds'])
                pval = hypergeom.sf(x - 1, len(pmn_M), n, N)
                if org not in org2enriched:
                    org2enriched[org]={}
                org2enriched[org][pt] = {}
                org2enriched[org][pt]['pval']=pval
                expected = (n * N) / len(pmn_M)
                fold = x / N
                org2enriched[org][pt]['coverage'] = fold
                org2enriched[org][pt]['total in pathway'] = N
                org2enriched[org][pt]['no. sampled'] = x
                org2enriched[org][pt]['expected'] = expected
                org2enriched[org][pt]['set']='PlantCyc'
                org2enriched[org][pt]['cat'] = pathway2cat[pt]
                org2enriched[org][pt]['cls'] = path2class[pt]

                if pval<0.05:
                    if org not in org2numenriched:
                        org2numenriched[org]=set()
                    org2numenriched[org].add(pt)
                    '''if pt not in ep2org:
                        ep2org[pt]=set()
                    ep2org[pt].add(org)'''
            else:
                bad.add(pt)
# convert results to dataframe and plot
eptcount=pd.DataFrame()

# create an org to pathway summary

orgpathnums=pd.DataFrame()
for org in org2pathways:
    kprim=0
    pprim=0
    ksec=0
    psec=0
    if 'kegg' in org2pathways[org]:
        kpaths=len(org2pathways[org]['kegg'])
        # calculate the primary and secondary coverage
        for path in org2pathways[org]['kegg']:
            if 'Prodigiosin biosynthesis' in path:
                continue
            if path in path2class:
                if 'Primary' in path2class[path]:
                    kprim+=1
                if 'Secondary' in path2class[path]:
                    ksec+=1

    else:
        kpaths=0
    if 'pmn' in org2pathways[org]:
        pmnpaths=len(org2pathways[org]['pmn'])
        for path in org2pathways[org]['pmn']:
            if path in path2class:
                if 'Primary' in path2class[path]:
                    pprim+=1
                if 'Secondary' in path2class[path]:
                    psec+=1
    else:
        pmnpaths=0
    totalp = kpaths+pmnpaths
    values_to_add = {'KEGG': kpaths, 'PlantCyc': pmnpaths, 'total': totalp,'KEGG primary':kprim,'KEGG secondary': ksec,
                     'PlantCyc primary': pprim,'PlantCyc secondary':psec}
    row_to_add = pd.Series(values_to_add, name=org)
    orgpathnums = orgpathnums.append(row_to_add)

orgpathnums.to_csv('gen_org2path.csv')


for org in org2numenriched:

    values_to_add = {'No. enriched pathways per plant':len(org2numenriched[org])}
    row_to_add = pd.Series(values_to_add, name=org)
    eptcount = eptcount.append(row_to_add)
    for path in totpaths:
        # fill out the organism dictionary in preparation for the plotting
        if path not in org2enriched[org] and (path in genPathPMN or path in genPathKEGG):
            org2enriched[org][path]={}
            org2enriched[org][path]['coverage']=0
            if path in genPathKEGG:
                org2enriched[org][path]['total in pathway'] = len(genPathKEGG[path]['compounds'])
                org2enriched[org][path]['set']='KEGG'
                org2enriched[org][path]['cat'] = pathway2cat[path]
                org2enriched[org][path]['cls'] = path2class[path]
            if path in genPathPMN:
                org2enriched[org][path]['total in pathway'] = len(genPathPMN[path]['compounds'])
                org2enriched[org][path]['set']='PlantCyc'
                org2enriched[org][path]['cat'] = pathway2cat[path]
                org2enriched[org][path]['cls'] = path2class[path]
            org2enriched[org][path]['no. sampled'] = 0

rows=[]
for org in org2enriched:
    data_row=org2enriched[org]
    for row in data_row:
        data_row[row]['plant']=org
        data_row[row]['pathway']=row

        rows.append(data_row[row])

orgsdf=pd.DataFrame(rows)

#try bonferoni
p_adjusted = multipletests(orgsdf['pval'], method='bonferroni')
orgsdf['bonferroni']=p_adjusted[1]

# get the primary vs. secondary for the general plants
tcounts=orgsdf.groupby('cls').count()
tpathcounts=orgsdf['pathway'].unique()
tprim=orgsdf[orgsdf['cls']=='Primary']
tsec=orgsdf[orgsdf['cls']=='Secondary']
tprimu=tprim['pathway'].unique()
tsecu=tsec['pathway'].unique()
prim=len(tprimu)/len(tpathcounts)
seco=len(tsecu)/len(tpathcounts)

sorgdf=orgsdf[orgsdf['bonferroni']<=0.05]
noten=orgsdf[orgsdf['bonferroni']>0.05]
# need to count not enriched per org and not enriched 2 orgs
stpathcounts=sorgdf['pathway'].unique()
stprim=sorgdf[sorgdf['cls']=='Primary']
stsec=sorgdf[sorgdf['cls']=='Secondary']
stprimu=stprim['pathway'].unique()
stsecu=stsec['pathway'].unique()
sprim=len(stprimu)/len(stpathcounts)
sseco=len(stsecu)/len(stpathcounts)


# get the corn dataframe and count unique pathways, also divide to primary and secondary
corn=orgsdf[orgsdf['plant']=='Zea mays']
cpaths=corn['pathway'].unique()
ccounts=corn.groupby('cls').count()
tprim=corn[corn['cls']=='Primary']
tsec=corn[corn['cls']=='Secondary']
tprimu=tprim['pathway'].unique()
tsecu=tsec['pathway'].unique()
prim=len(tprimu)/len(tpathcounts)
seco=len(tsecu)/len(tpathcounts)



enpaths=corn[corn['bonferroni']<=0.05]
notenpaths=corn[corn['bonferroni']>0.05]
encounts=enpaths.groupby('cls').count()

stpathcounts=enpaths['pathway'].unique()
stprim=enpaths[enpaths['cls']=='Primary']
stsec=enpaths[enpaths['cls']=='Secondary']
stprimu=stprim['pathway'].unique()
stsecu=stsec['pathway'].unique()
sprim=len(stprimu)/len(stpathcounts)
sseco=len(stsecu)/len(stpathcounts)



plt.figure(figsize=(25, 10))
sns.set(font_scale=1.5)
p=sns.boxplot(x="cat", y="bonferroni", data=corn)
p.set_xticklabels(p.get_xticklabels(), rotation=90, fontsize=16)
plt.axhline(y=0.05, color='black', linestyle='--')
p.set(xlabel='Category', ylabel='P-value (Bonferroni)')
#plt.xticks(catlist)
plt.tight_layout()
plt.show()



# add a boxplot of the enriched categories
catlist=list(orgsdf['cat'].unique())
sorgdf=orgsdf[orgsdf['bonferroni']<=0.05]
plt.figure(figsize=(25, 10))
sns.set(font_scale=1.5)
p=sns.boxplot(x="cat", y="bonferroni", data=orgsdf)
p.set_xticklabels(p.get_xticklabels(), rotation=90, fontsize=16)
plt.axhline(y=0.05, color='black', linestyle='--')
p.set(xlabel='Category', ylabel='P-value (Bonferroni)')
#plt.xticks(catlist)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 10))
p=sns.boxplot(x="cls", y="bonferroni", data=orgsdf)
p.set_xticklabels(p.get_xticklabels(), rotation=90, fontsize=16)
plt.axhline(y=0.05, color='black', linestyle='--')
p.set(xlabel='Metabolism', ylabel='P-value (Bonferroni)')
#plt.xticks(catlist)
plt.tight_layout()
plt.show()


# create an updated ep2org
ep2org={}
for org in org2numenriched:
    # update the number for the adjusted pvalue
    odf=orgsdf[orgsdf['plant']==org]
    odf=odf[odf['bonferroni']<0.05]
    pathways=odf['pathway'].unique()
    for p in pathways:
        if p not in ep2org:
            ep2org[p]=set()
        ep2org[p].add(org)



orgsdf.to_csv('org_pathway_enrichment_enriched_corr.csv')
orgsdf=orgsdf.dropna()
# get the corrected lower than 0.05 set
enbondf=orgsdf[orgsdf['bonferroni']<0.05]
# try out just kegg
orgsdf_pmn=enbondf[enbondf['set']=='PlantCyc']
orgsdf_kegg=enbondf[enbondf['set']=='KEGG']

#plotPCA(orgsdf_pmn,'orgsdf_pmn')
'''plotPCA(orgsdf_kegg,'orgsdf_kegg')
plotPCA(enbondf,'total')'''







corndf=enbondf[enbondf['plant']=='Zea mays']

corndf=corndf[corndf['coverage']>0]
colors = ['#BF4D5C',"#FF0B04", "#4374B3"]
# plot corns coverage distribution
#distribution
plt.figure(figsize=(18, 10))
sns.set(font_scale=2.5)
g=sns.distplot(corndf['coverage'],kde_kws=dict(linewidth=5), color='#BF4D5C')
plt.tight_layout()
plt.show()

corndf_pmn=corndf[corndf['set']=='PlantCyc']
corndf_kegg=corndf[corndf['set']=='KEGG']
#plotPCA(corndf_pmn,'corndf_pmn')
#plotPCA(corndf_kegg,'corndf_kegg')
# plot a pca of the corn subset - once for kegg and once for plantcyc
# KEGG








bon2org={}


ptorgcount=pd.DataFrame()
singles=pd.DataFrame()
for pt in ep2org:
    values_to_add = {'No. plants per enriched pathway': len(ep2org[pt])}
    row_to_add = pd.Series(values_to_add, name=pt)
    ptorgcount = ptorgcount.append(row_to_add)
    # single out the enriched for a single organism
    if len(ep2org[pt])==1:
        org=next(iter(ep2org[pt]))
        if org in org2enriched and pt in org2enriched[org] and org2enriched[org][pt]['pval']>=0.05:
            continue
        values_to_add = {'plant': org, 'Pathway':pt,'Coverage': org2enriched[org][pt]['coverage'],'pval':org2enriched[org][pt]['pval']}
        row_to_add = pd.Series(values_to_add, name=org)
        singles = singles.append(row_to_add)
singles.to_csv('singles.csv')
# plot distribution
#distribution
plt.figure(figsize=(18, 6))
sns.set(font_scale=2.5)
g=sns.distplot(eptcount['No. enriched pathways per plant'], color='dodgerblue',kde_kws=dict(linewidth=5))
plt.tight_layout()
plt.show()
#distribution
plt.figure(figsize=(18, 6))
sns.set(font_scale=2.5)
g=sns.distplot(ptorgcount['No. plants per enriched pathway'], color='forestgreen',kde_kws=dict(linewidth=5))
plt.tight_layout()
plt.show()

# plot the singles in a bar plot
cp1=['#ffd429','#73AF59','#fbec5d','#ede6d9','wheat','palegoldenrod','purple','burlywood','seagreen','olivedrab',
     '#daa520','#ffb6c1','#008000','orange','rosybrown','gold','greenyellow','red','brown','#ffa500',
     'darkorange','y','yellowgreen','crimson','ivory','mediumvioletred','firebrick','darkgreen','chocolate',
     'tomato','darkseagreen','khaki']
plt.figure(figsize=(30, 12))
sns.set_palette(sns.color_palette(cp1))
sns.set(font_scale=2.5)
p=sns.barplot(x='Pathway',y='Coverage', data=singles, hue=singles.index,dodge=False,palette=sns.color_palette(cp1),edgecolor='gray')
p.set_xticklabels(p.get_xticklabels(), rotation=90, fontsize=10)
#p.set(xticks=[])

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=5, mode="expand", borderaxespad=0.)
plt.tight_layout()
plt.savefig("f2c.svg")
plt.show()

# make a filtered version with the seven plants I'm talking about in the manuscript
orgshlist=['Vitis vinifera','Hordeum vulgare subsp. Vulgare','Ricinus communis','Capsicum annuum',
           'Solanum lycopersicum','Sorghum bicolor','Zea mays']
singlesf=pd.DataFrame()
for org in orgshlist:
    tdf=singles[singles['plant']==org]
    singlesf=singlesf.append(tdf)
cp2=['#6f2da8','#84563c','#cb4154','red','tomato','palegoldenrod','#fbec5d','burlywood','seagreen','olivedrab',
     '#daa520','#ffb6c1','#008000','orange','rosybrown','gold','greenyellow','red','brown','#ffa500',
     'darkorange','y','yellowgreen','crimson','ivory','mediumvioletred','firebrick','darkgreen','chocolate',
     'tomato','darkseagreen','khaki']
plt.figure(figsize=(30, 15))
sns.set_palette(sns.color_palette(cp1))
sns.set(font_scale=2.5)
p=sns.barplot(x='Pathway',y='Coverage', data=singlesf, hue=singlesf.index,dodge=False,palette=sns.color_palette(cp2),edgecolor='gray')
p.set_xticklabels(p.get_xticklabels(), rotation=90, fontsize=26)
#p.set(xticks=[])

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0.)
plt.tight_layout()
plt.savefig("f2cf.svg")
plt.show()

pickle.dump( org2pathways, open( "org2pathways.pkl", "wb" ) )



# try a different clustering on the pathways - hirarchical
#create the mtx
pathlist=list(org2enriched['Malus domestica'])
orglist=orgsdf.plant.unique()
mtx=pd.DataFrame()

for org in org2enriched:
    tlist=[]
    for pt in pathlist:
        if pt in org2enriched[org]:
            if 'pval' in org2enriched[org][pt]:
                tlist.append(org2enriched[org][pt]['pval'])
            else:
                tlist.append(10)

    mtx[org]=tlist


# I have the matrix - let's cluster
import scipy.cluster.hierarchy as shc
plt.figure(figsize=(10, 7))
dend = shc.dendrogram(shc.linkage(mtx, method='ward'))
plt.axhline(y=300, color='black', linestyle='--')
plt.show()



from sklearn.cluster import AgglomerativeClustering
cluster = AgglomerativeClustering(n_clusters=4, affinity='euclidean', linkage='ward')
cluster.fit_predict(mtx)

# try k-means
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

kmeans = KMeans( init="random",   n_clusters=4, n_init=10,max_iter=300,random_state=42 )
kmeans.fit(mtx)
'''plt.figure(figsize=(10, 7))
plt.scatter(mtx, c=cluster.labels_)
plt.show()'''
labels=cluster.labels_
mtx['cluster']=labels

mtx['kmeans_cluster']=kmeans.labels_
resdf=pd.DataFrame()
resdf['agglo_cluster']=labels
resdf['kmeans_cluster']=kmeans.labels_
resdf['pathway']=pathlist
resdf.to_csv('clustering results_paths')

# make a 3d umap plot
from umap import UMAP
umap_3d = UMAP(n_components=3, init='random', random_state=0)
proj_3d = umap_3d.fit_transform(mtx)
fig_3d = px.scatter_3d(
    proj_3d, x=0, y=1, z=2,
    color=mtx.cluster, labels={'color': 'cluster'}
)
fig_3d.update_traces(marker_size=5)
fig_3d.show()
print ('Done')