''' creating figure 2 for the genomics papaer about corn'''


# imports
import math
import re

import venn
from bioservices import *
import json
from scipy.stats import hypergeom
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import pickle
from tqdm import tqdm
# --------------- functions ----------------

def cleanhtml(raw_html):
    # html tag cleaner copied from stack overflow, cleans non tag entities as well
  cleanr = re.compile('<.*?>|&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-f]{1,6});')
  cleantext = re.sub(cleanr, '', raw_html)
  return cleantext

def pathClassify(pathway):
    '''classifies the compounds in a pathway according to their sources if possible - fdb,dfc, usda or genomics'''

    if pathway not in genPath:
        pathway = pathway.replace('"', '')
        pathway=pathway.replace('plants, ','')
        pathway=cleanhtml(pathway)
        '''if any(x in pathway for x in ['Biosynthesis of secondary metabolites','Metabolic pathways',
                                      'Microbial metabolism in diverse environments','palmitate biosynthesis (type II fatty acid synthase)',
                                      'octanoyl-[acyl-carrier protein] biosynthesis (mitochondria, yeast)',
                                      'phospholipid remodeling (phosphatidylcholine, yeast)',
                                      'acetyl-CoA biosynthesis from citrate','pyruvate decarboxylation to acetyl CoA I',
                                      'mevalonate pathway I (eukaryotes and bacteria)','fatty acid biosynthesis initiation (type II)'
                                      ,'fatty acid beta-oxidation II (plant peroxisome)','L-Ndelta-acetylornithine biosynthesis',
                                      'phospholipid remodeling (phosphatidate, yeast)','beta-alanine biosynthesis II',
                                      'Fatty acid metabolism','Carbon metabolism','benzoate biosynthesis III (CoA-dependent, non-beta-oxidative)',
                                      'Degradation of aromatic compounds','Biosynthesis of amino acids','2-Oxocarboxylic acid metabolism',
                                      'lipid IVA biosynthesis (E. coli)','reactive oxygen species degradation','NAD salvage (plants)',
                                      'diphthamide biosynthesis I (archaea)','folate transformations II (plants)',
                                      'S-methyl-5-thio-alpha-D-ribose 1-phosphate degradation I','L-glutamate biosynthesis III',
                                      'beta-alanine biosynthesis I',
                                      'UDP-alpha-D-glucuronate biosynthesis (from myo-inositol)','alpha-carotene biosynthesis',
                                      'molybdopterin biosynthesis','ethene biosynthesis I (plants)','Kdo transfer to lipid IVA I (E. coli)',
                                      'L-proline degradation I','L-proline biosynthesis I (from L-glutamate)',
                                      'UDP-alpha-D-galacturonate biosynthesis II (from D-galacturonate)']):'''
        return

    for comp in genPath[pathway]['compounds']:
        # get the inchi key using the name of the compound
        if comp in name2ik:
            ik=name2ik[name].strip()
            if ik in corntot:
                for src in mdm[ik]['source']:
                    if 'NDM' in src:
                        continue
                    if pathway not in gencornpathways:
                        gencornpathways[pathway]={}
                    if comp not in gencornpathways[pathway]:
                        gencornpathways[pathway][comp]=set()
                    gencornpathways[pathway][comp].add(src)
    return


def writetxt(itera, fname):
    outf=open(fname,'w')
    for item in itera:
        if isinstance(item,float):
            pass
        if  item in ik2name and isinstance(ik2name[item],float):
            pass
        if item in ik2name and not isinstance(item,float):
            outf.write(str(item)+'\t'+str(ik2name[item])+'\n')
        else:
            outf.write(str(item) + '\n')
    outf.close()


# ---------------- main ----------------------
mdm=pickle.load(open('mdm_040621.pkl','rb'))
ndmdf=pd.read_csv('NDM_MasterJan21.csv')
ndmdf.set_index('InChIKey')
ndm=ndmdf.set_index('InChIKey').T.to_dict('list')
#ik2name={}

#dump it
#pickle.dump( ik2name, open( "ik2name_jun21.pkl", "wb" ) )

# load it
ik2name=pickle.load(open('ik2name_jun21.pkl','rb'))


# first losd the corn numbers

orgs=pickle.load(open('orgkeys_new.pkl','rb'))
# create the stacked bar of all organisms
tfdb=set()
tusda=set()
tdfc=set()
totkegg=set()
totpmn=set()
totothers=set()
totonlyexpri=set()
df=pd.DataFrame()
# for the overlaps I'm going to use the first block of the inchi to align with performance mcalculations
# that means I need to generate a fb sets for all
for org in orgs:
    total=0
    tot=set()
    genomics=set()
    if 'kegg' in orgs[org]:
        genomics=genomics|orgs[org]['kegg']
        kegg=len(orgs[org]['kegg'])
        tot=tot|genomics
        totkegg=totkegg|orgs[org]['kegg']
    else:
        kegg=0
    if 'pmn' in orgs[org]:
        genomics=genomics|orgs[org]['pmn']
        pmn=len(orgs[org]['pmn'])
        tot=tot|genomics
        totpmn=totpmn|orgs[org]['pmn']
    else:
        pmn=0
    gens=len(genomics)
    if 'foodb' in orgs[org]:
        fdb=len(orgs[org]['foodb'])
        tfdb=tfdb|orgs[org]['foodb']
        tot=tot|orgs[org]['foodb']
        totothers=totothers|orgs[org]['foodb']
    else:
        fdb=0
    if 'usda' in orgs[org]:
        usda=len(orgs[org]['usda'])
        tusda=tusda|orgs[org]['usda']
        tot=tot|orgs[org]['usda']
        totothers=totothers|orgs[org]['usda']
    else:
        usda=0
    if 'dfc' in orgs[org]:
        dfc=len(orgs[org]['dfc'])
        tdfc=tdfc|orgs[org]['dfc']
        tot=tot|orgs[org]['dfc']
        totothers=totothers|orgs[org]['dfc']
    else:
        dfc=0
    if 'experiment' in orgs[org]:
        experiments=len(orgs[org]['experiment'])
        tot=tot|orgs[org]['experiment']
        for atr in orgs[org]:
            if atr in orgs[org]:
                totonlyexpri = totonlyexpri | orgs[org][atr]

    else:
        experiments=0
    total=fdb+dfc+experiments+kegg+pmn
    values_to_add = {'1_PlantCyc': pmn, '2_KEGG':kegg, '3_FooDB': fdb, '5_USDA': usda,
                     '4_DFC': dfc, '6_Experiments': experiments,'Plant': org,'Total': total}
    row_to_add = pd.Series(values_to_add, name=org)
    df = df.append(row_to_add)
df=df.sort_values(by='Total')
df=df.drop(columns=['Total'])
# plot
# plot a Stacked Bar Chart using matplotlib
print(df.mean())

print(df.std())
print(df.median())
df.plot(
  x = 'Plant',
  kind = 'barh',
  stacked = True,
  figsize=(40,30),
  mark_right = True,
fontsize=30,
color=['forestgreen','yellowgreen','gold','red','dodgerblue','orange'])
plt.ylabel('Plant', fontsize=40, fontweight='bold')
plt.xlabel('No. Compounds', fontsize=40,fontweight='bold')
plt.legend(fontsize=30,markerscale=2.)
plt.tight_layout()
plt.show()


plt.figure(figsize=(12,12))
sns.set(font_scale=2)
out=venn3([totothers, totkegg,totpmn], set_labels=('FooDB, DFC and USDA', 'KEGG','PMN'),
      set_colors=('grey', 'yellowgreen','forestgreen'), alpha=0.7)
plt.tight_layout()
plt.show()


# get the sets for all the sources
outdf=pd.DataFrame()
corng=orgs['zea mays']['kegg']|orgs['zea mays']['pmn']
writetxt(corng,'corn_genomics_subset.txt')
cornfdb=orgs['zea mays']['foodb']
writetxt(cornfdb,'corn_fdb_subset.txt')

corndfc=orgs['zea mays']['dfc']
writetxt(corndfc,'corn_dfc_subset.txt')
cornusda=orgs['zea mays']['usda']
writetxt(cornusda,'corn_usda_subset.txt')
cornexp=orgs['zea mays']['experiment']
writetxt(cornexp,'corn_exp_subset.txt')

# add corn kinetic keys
cornk=pd.read_csv('kcorn.csv')
cornkd=cornk.set_index('inchi key').to_dict('index')
cornkinteics=set(list(cornkd.keys()))
#manual additions from curation
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
corng=corng|set(cur)
corncheck=cornfdb|corndfc|cornusda|corng|cornexp

for ik in corncheck:
    if ik in corng:
        g=1
    else:
        g=0
    if ik in cornfdb:
        fd=1
    else:
        fd=0
    if ik in corndfc:
        df=1
    else:
        df=0
    if ik in cornusda:
        us=1
    else:
        us=0
    if ik in cornexp:
        cexp=1
    else:
        cexp=0
    # update the dataframe
    values_to_add = {'genomics': g, 'foodb': fd, 'dfc': df, 'usda': us, 'experiments': cexp}
    row_to_add = pd.Series(values_to_add, name=ik)
    outdf = outdf.append(row_to_add)

outdf.to_csv('corn_subsets_binary.csv')
plt.figure(figsize=(12,12))
sns.set(font_scale=2.5)
out=venn3([cornfdb, corndfc,cornusda], set_labels=('FooDB', 'DFC','USDA'),
      set_colors=('yellowgreen', 'dodgerblue','tomato'), alpha=0.7)
plt.show()


# plot the overlap of the total sets
labels = venn.get_labels([corng,cornusda,cornfdb,corndfc, cornexp], fill=['number'])
fig, ax = venn.venn5(labels, names=['Genomes', 'USDA', 'FooDB', 'DFC','Experiments'])
fig.show()

out=venn3([corng,cornkinteics,cornexp], set_labels=('Genomics','Kinetics','Experiments'),
      set_colors=('yellowgreen','dodgerblue','orange'), alpha=0.7)
plt.show()


# get the first block version of orgkeys
orgsfb=pickle.load(open('orgkeys_fb.pkl','rb'))
# load the kintics to add to significance analysis
orgdgscores=pickle.load(open('kineticScores.pkl','rb'))

orgdgscoresnew={}
for org in orgdgscores:
    orgl=org.lower()
    orgdgscoresnew[orgl]=orgdgscores[org]

orgdgscores=orgdgscoresnew

# attach names to the keys in the kinetics analysis
for org in orgdgscores:
    if 'experiment' in orgs[org]:
        f1=open(org+'_kinetic_scores.txt','w')
        f2=open(org+'avg_kinetic_scores.txt','w')
        kik2score={}
        kik2avgscore={}
        for fb in orgdgscores[org].keys():
            for ik in orgdgscores[org][fb]['inchikey']:
                if ik in orgs[org]['experiment']:
                    exp=1
                else:
                    exp=0
                kik2score[ik]=orgdgscores[org][fb]['sum']
                if ik in ik2name:
                    f1.write(ik+'\t'+str(ik2name[ik])+'\t'+str(orgdgscores[org][fb]['sum'])+'\t'+str(exp)+'\n')
                    f2.write(ik + '\t' + str(ik2name[ik]) + '\t' + str(orgdgscores[org][fb]['average'])+'\t' +str(exp)+ '\n')
                else:
                    f1.write(ik + '\t\t' +str(orgdgscores[org][fb]['sum'])+'\t'+ str(exp)+ '\n')
                    f2.write(ik + '\t\t' + str(orgdgscores[org][fb]['average']) +'\t'+str(exp)+ '\n')
f1.close()
f2.close()


totonlyexpri=set()
for org in orgsfb:
    if 'experiment' in orgsfb[org]:
        for atr in orgsfb[org]:
            if atr in orgsfb[org]:
                totonlyexpri = totonlyexpri | orgsfb[org][atr]



# calculate the uniqur genome vs. experiments and existing for all plants, then plot
ovdf=pd.DataFrame()
for org in orgs:
    if 'experiment' in orgs[org]:
        if 'kegg' not in orgs[org]:
            corng = orgs[org]['pmn']
        elif 'pmn' not in orgs[org]:
            corng = orgs[org]['kegg']
        else:
            corng = orgs[org]['kegg'] | orgs[org]['pmn']
        #manually curated inchi key discrepencies fix as a result of keys coming from different sources
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
        fbcur=set()
        for c in cur:
            cl=c.split('-')
            fbcur.add(cl[0])
        corng = corng | set(cur)
        cornfdb = orgs[org]['foodb']
        corndfc = orgs[org]['dfc']
        cornusda = orgs[org]['usda']
        cornexp = orgs[org]['experiment']
        # first block versions of these
        if 'kegg' not in orgsfb[org]:
            corngfb = orgsfb[org]['pmn']
        elif 'pmn' not in orgsfb[org]:
            corngfb = orgsfb[org]['kegg']
        else:
            corngfb = orgsfb[org]['kegg'] | orgsfb[org]['pmn']
        corngfb = corngfb | set(fbcur)
        cornfdbfb=orgsfb[org]['foodb']
        corndfcfb=orgsfb[org]['dfc']
        cornusdafb=orgsfb[org]['usda']
        cornexpfb=orgsfb[org]['experiment']
        expandgenunique = (corng - cornfdb - cornusda - corndfc) & cornexp
        totexisting = set()
        eonlyusda = (cornusda - cornfdb - corndfc - corng) & cornexp
        totexisting = totexisting | eonlyusda

        eusdafdb = (cornusda - corng - corndfc) & cornexp & cornfdb
        totexisting = totexisting | eusdafdb

        eonlyfdb = (cornfdb - corng - cornusda - corndfc) & cornexp
        totexisting = totexisting | eonlyfdb

        edfcfdb = (corndfc - corng - cornusda) & cornexp & cornfdb
        totexisting = totexisting | edfcfdb

        eonlydfc = (corndfc - corng - cornusda - cornfdb) & cornexp
        totexisting = totexisting | eonlydfc

        # formulate the hypergeometric test
        # all we know about the organisms that have experiments
        M = len(totonlyexpri)
        # the existing set
        eset=cornfdbfb|corndfcfb|cornusdafb
        N = len(eset)
        # the genomics set
        Nc = len(corngfb)
        # the positive kinetics set
        # gather the positive scored fb
        orgpos = set()
        totkinetics=set()
        orgikpos = set()
        for fb in orgdgscores[org]:
            totkinetics.add(fb)
            if orgdgscores[org][fb]['sum'] > 0:
                orgpos.add(fb)
                for ik in orgdgscores[org][fb]['inchikey']:
                    orgikpos.add(ik)
        Nkt=len(orgpos)
        Nktt=len(totkinetics)
        # experiments
        n = len(orgsfb[org]['experiment'])
        # overlap between existing and experiments
        overlap_e =  eset & orgsfb[org]['experiment']
        x = len(eset & orgsfb[org]['experiment'])

        # overlap between the experiments and the genomics organism set
        toverlap = orgsfb[org]['experiment'] & corngfb
        xc = len(orgsfb[org]['experiment'] & corngfb)

        # overlap between the experiments and the positive kinetics organism set
        tkoverlap = orgsfb[org]['experiment'] & orgpos
        xct = len(orgsfb[org]['experiment'] & orgpos)

        # overlap between the experiments and the total kinetics organism set
        ttkoverlap = orgsfb[org]['experiment'] & totkinetics
        xctt = len(orgsfb[org]['experiment'] & totkinetics)

        # test existing
        pval = hypergeom.sf(x - 1, M, n, N)
        expected = (n * N) / M
        fold = x / N



        # test genomics
        pvalc = hypergeom.sf(xc - 1, M, n, Nc)
        expectedc = (n * Nc) / M
        fold = xc / Nc

        # test  positive kinetics
        pvalk = hypergeom.sf(xct - 1, M, n, Nkt)
        expectedk = (n * Nkt) / M
        fold = xct / Nkt

        # test  total kinetics
        pvaltk = hypergeom.sf(xctt - 1, M, n, Nktt)
        expectedtk = (n * Nktt) / M
        fold = xctt / Nktt

        # update the dataframe
        values_to_add = {'plant': org,'unique genomics': len(expandgenunique),'unique existing knowledge':len(totexisting),
                         'Overlap existing': x, 'P-value existing': pval, 'Expected overlap existing': expected,
                         'Overlap genomics': xc,
                         'Expected overlap genomics': expectedc, 'P-value genomics': pvalc,
                         'Expected overlap kinetics': expectedtk, 'Overlap kinetics':xctt, 'P-value kinetics': pvaltk
                         }
        row_to_add = pd.Series(values_to_add, name=org)
        ovdf = ovdf.append(row_to_add)


print(ovdf.mean())
print(ovdf.std())
ovdf=ovdf.sort_values(by=['unique genomics'])
# plot


color={'unique genomics':'forestgreen',
       'unique existing knowledge': 'grey'}
ovdf.plot(x="plant", y=["unique genomics", "unique existing knowledge"], kind="bar", color=color, figsize=(20,15))
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.ylabel("No. of compounds")
plt.tight_layout()
plt.show()


sns.set(font_scale=3)
color={'Overlap existing':'steelblue',
       'Expected overlap existing': 'grey',
        'Overlap genomics':'yellowgreen',
        'Expected overlap genomics': 'darkseagreen',
       'Expected overlap kinetics': 'orange',
        'Overlap kinetics':'goldenrod'
       }
ovdf.plot(x="plant", y=["Expected overlap existing","Overlap existing",'Expected overlap genomics','Overlap genomics',
                        'Expected overlap kinetics','Overlap kinetics'], kind="bar", color=color,width=0.8, figsize=(30,20))
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.ylabel("No. of compounds")
plt.tight_layout()
plt.show()

sns.set(font_scale=2.5)
plt.figure(figsize=(15, 15))
f1=sns.scatterplot(data=ovdf, x="plant", y="P-value existing", color='steelblue', s=100)
f2=sns.scatterplot(data=ovdf, x="plant", y="P-value genomics", color='yellowgreen', s=100)
f2=sns.scatterplot(data=ovdf, x="plant", y="P-value kinetics", color='goldenrod',s=100)
f1.set(yscale="log")
plt.axhline(y=0.05, color='black', linestyle='--')
plt.draw()
f1.set_xticklabels(f1.get_xticklabels(), rotation=90)

plt.tight_layout()
plt.show()

ik2type=pickle.load(open('ik2prim_sec.pkl','rb'))
# get the info on the keys
expudf=pd.DataFrame()
for ik in onlyexp:
    ikser=ndmdf[ndmdf['InChIKey']==ik]
    if ik in ik2type:
        ikser['type']=ik2type[ik]
    else:
        ikser['type']='unknown'
    expudf=expudf.append(ikser)

expudf.to_csv('only_exp_corn.csv')

expudf=pd.DataFrame()
for ik in expandgenunique:
    ikser=ndmdf[ndmdf['InChIKey']==ik]
    if ik in ik2type:
        if len(ik2type[ik])==2:
            ikser['type'] = 'secondary'
        else:
            ikser['type']='primary'
    else:
        ikser['type']='unknown'
    expudf=expudf.append(ikser)
expudf.to_csv('exp_gen_corn.csv')




global corntot
corntot=corng|cornfdb|corndfc|cornusda|cornexp
# create data
names = 'Genomics', 'FooDB', 'DFC', 'USDA', 'Experiments'
size = [len(corng), len(cornfdb), len(corndfc), len(cornusda), len(cornexp)]


# Create a circle for the center of the plot
sns.set(font_scale=1.6)
my_circle = plt.Circle((0, 0), 0.5, color='white')
# Give color names
plt.pie(size, labels=names, colors=['yellowgreen','gold','red','dodgerblue','orange'],
        wedgeprops = { 'linewidth' : 7, 'edgecolor' : 'white' }, autopct='%1.0f%%', pctdistance=1.2, labeldistance=1.4)
p=plt.gcf()
p.gca().add_artist(my_circle)
plt.tight_layout()
plt.show()

# now we need a pathway example that will show the coverage from different sources.For that I will use MDM
# load mdm

name2ik={}
ndmtots=pickle.load(open('NDMTotalsMArch1.pkl','rb'))
cornkinetics=pickle.load(open('kcorn.pkl','rb'))
for ik in mdm:
    # add fdb,dfc and usda sources to mdm
    for tot in ndmtots:
        if ik in ndmtots[tot]:
            mdm[ik]['source'].append(tot)
    if 'name' in mdm[ik]:
        name=mdm[ik]['name'].lower()
        if name not in name2ik:
            name2ik[name]=ik
    if 'kegg' in mdm[ik]:
        if isinstance(mdm[ik]['kegg'],list):
            name=mdm[ik]['kegg'][0]
        else:
            name=mdm[ik]['kegg']
        if name not in name2ik:
            name2ik[name]=ik
    if 'plantcyc' in mdm[ik]:
        name = mdm[ik]['plantcyc'].lower()
        if name not in name2ik:
            name2ik[name] = ik
# go over the general pathways and collect compounds that are not in the genomics set but in a pathway
global genPath
genPath=pickle.load(open('genpaths.pkl','rb'))
global gencornpathways
gencornpathways={}
cornpaths={}
pdf=pd.DataFrame()
for ik in tqdm(mdm):


    if 'orgs' in mdm[ik] and 'Zea mays' in mdm[ik]['orgs']:
        cornd=mdm[ik]['orgs']['Zea mays']
        for r in cornd:
            # get pathways
            pa = cornd[r]
            for p in pa:
                if not isinstance(pa[p], dict):
                    # pmn
                    if 'no mapping' in p:
                        continue
                    # classify full pathway
                    pathClassify(p)
                    # add pathway to list
                    if p not in cornpaths:
                        cornpaths[p]={}
                    # add compound name and sources
                    pmncomp=mdm[ik]['plantcyc'].lower()
                    if pmncomp not in cornpaths[p]:
                        cornpaths[p][pmncomp]=set()
                    for src in mdm[ik]['source']:
                        if 'NDM' in src:
                            continue
                        cornpaths[p][pmncomp].add(src)
                    # check for dfc and usda
                    if ik in corndfc:
                        cornpaths[p][pmncomp].add('dfc')
                        mdm[ik]['source'].append('dfc')
                    if ik in cornusda:
                        cornpaths[p][pmncomp].add('usda')
                        mdm[ik]['source'].append('usda')
                    if ik in cornfdb:
                        cornpaths[p][pmncomp].add('fdb')
                        mdm[ik]['source'].append('fdb')

                else:
                    # a kegg situation, r is infact an enzyme
                    for path in pa[p]:
                        if 'no mapping' in path:
                            continue
                        # classify general pathway
                        pathClassify(path)
                        # add pathway to list
                        if path not in cornpaths:
                            cornpaths[path] = {}
                        # add compound name and sources
                        if isinstance(mdm[ik]['kegg'], list):
                            keggcomp = mdm[ik]['kegg'][0].lower()
                        else:
                            keggcomp=mdm[ik]['kegg'].lower()
                        if keggcomp not in cornpaths[path]:
                            cornpaths[path][keggcomp] = set()
                        for src in mdm[ik]['source']:
                            if 'NDM' in src:
                                continue
                            cornpaths[path][keggcomp].add(src)
                        # check for dfc and usda
                        if ik in corndfc:
                            cornpaths[path][keggcomp].add('dfc')
                            mdm[ik]['source'].append('dfc')
                        if ik in cornusda:
                            cornpaths[path][keggcomp].add('usda')
                            mdm[ik]['source'].append('usda')
                        if ik in cornfdb:
                            cornpaths[path][keggcomp].add('fdb')
                            mdm[ik]['source'].append('fdb')



# convert to managable dataframe
for path in cornpaths:
    if path in gencornpathways:
        cornpaths.update(gencornpathways[path])
    for comp in cornpaths[path]:
        genomics,fdb,dfc,usda=0,0,0,0
        if 'KEGG' in cornpaths[path][comp]:
            genomics=1
        if 'PlantCyc' in cornpaths[path][comp]:
            genomics=1
        if 'fdb' in cornpaths[path][comp]:
            fdb=1
        if 'dfc' in cornpaths[path][comp]:
            dfc=1
        if 'usda' in cornpaths[path][comp]:
            usda=1

        values_to_add = {'pathway':path,'compound':comp,'genomics':genomics,'foodb':fdb,
                         'dfc':dfc,'usda': usda}
        row_to_add = pd.Series(values_to_add,name=path)
        pdf = pdf.append(row_to_add)


'''# plot heatmap of each pathway and choose the most representing
for path in tqdm(cornpaths):
    psdf = pdf[pdf["pathway"] ==path]
    if 1.0 in psdf['usda'].values:
        if 'biosynthesis' in path:
            # plot heatmap
            plt.figure(figsize=(15, 10))
            sns.set(font_scale=1)
            sns.heatmap(psdf[['dfc', 'foodb', 'genomics', 'usda']], linewidths=.5, cmap="Greys")
            plt.tight_layout()
            #plt.show()'''



print('Done')