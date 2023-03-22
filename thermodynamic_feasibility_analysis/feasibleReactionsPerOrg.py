''' Add thermodynamic information to the reactions in MDM.
 Then try and identify what is an intermediate and what is an end product by using deltaG value for the available reactions
 finally try and compare to the experimental results.
 Here the calculation is done for a specific organism

 Add constraints to account for mass spec limitations:
 1. remove all compounds with wieght lower than 50.
 2. Use only the first block of the inchi key to disregard streochemistry. Masspec doesn't tell between stereoisomers.
 3. bug fixes - reaction counter and sum is wrong.

 Written by Shany Ofaim , the Barabasi lab CCNR Northeastern University 2021'''

# imports
import json
import pickle
import statistics
from time import sleep
from tqdm import tqdm
from bioservices import *
import pubchempy as pcp
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import re


# ------ functions ----------
def verifyIk(iklistsource, iklisttarget):
    # verify that the inchi key is correct and update the relevant dictionaries.
    # first check all sources for matches and then check with pubchem or chebi.
    fout = open('badiklist.txt', 'w')
    badik = 0
    for cnt, ik in enumerate(iklisttarget):
        print(str(cnt))
        if ik not in iklistsource:
            # inchi keys do not match in this case I'm ruling in favor of mdm. update relevant dictionaries
            badik += 1

            if ik in ik2kegg and ik2kegg[ik] in kegg2mdm:
                # found a match in mdm now replace the inchi key in the target list
                kname = ik2kegg[ik]
                # update the dictionary
                kegg2ik[kname] = kegg2mdm[kname]
            elif ik in ik2metacyc and ik2metacyc[ik] in meta2mdm:
                pname = ik2metacyc[ik]
                metacyc2ik[pname] = meta2mdm[pname]
            else:
                # try getting info from pubchem
                while True:
                    if ik in ik2name:
                        for na in ik2name[ik]:
                            try:
                                revcomp = pcp.get_compounds(na, 'name')
                                if revcomp:
                                    pcik = revcomp[0].inchikey
                                    if pcik not in ik:
                                        # check if ik gets a hit in pubchem
                                        testcomp = pcp.get_compounds(ik, 'inchikey')
                                        if not testcomp:
                                            fout.write(ik2cpd[ik] + '\t' + ik + '\t' + pcik + '\n')
                                    if pcik in mdm:
                                        if 'kegg' in mdm[pcik]:
                                            kegg2ik[mdm[pcik]['kegg']] = pcik
                                        else:
                                            if ik in ik2kegg:
                                                mdm[pcik]['kegg'] = ik2kegg[ik]
                                        if 'plantcyc' in mdm[pcik]:
                                            metacyc2ik[mdm[pcik]['plantcyc']] = pcik
                                        else:
                                            if ik in ik2metacyc:
                                                mdm[pcik]['plantcyc'] = ik2metacyc[ik]
                                        break

                                    else:
                                        notinmdm[ik] = na
                                        break
                                    break

                            except Exception as ex:
                                # If an error is raised it will print it while the code runs.
                                exc_type, exc_obj, exc_tb = sys.exc_info()
                                print(exc_type, exc_tb.tb_lineno)
                                print('Retrying in 0.5 seconds. Exception in connecting')
                                sleep(0.5)
                        break
                    else:
                        break
    fout.close()
    return badik


def getSubsProdsKEGG(equation, direction,org):
    # gets the substrates and products of a reaction from the equation and direction
    # reversible reactions will go from left to right
    if '<=>' in equation:
        eq = equation.split(' <=> ')
    elif '=>' in equation:
        eq = equation.split(' => ')
    else:
        eq = equation.split(' <= ')
    left = eq[0].split(' + ')
    right = eq[1].split(' + ')
    if '>' or '=' in direction:
        orgreacs[org][p]['substrates'] = left
        orgreacs[org][p]['products'] = right
    if '<' in direction:
        orgreacs[org][p]['substrates'] = right
        orgreacs[org][p]['products'] = left
    return


def getSubsProdsSEED(equation, direction,org):
    # gets the substrates and products of a reaction from the equation and direction
    # reversible reactions will go from left to right
    # add check on the inchi key from seed
    if '<=>' in equation:
        eq = equation.split(' <=> ')
    elif '=>' in equation:
        eq = equation.split(' => ')
    else:
        eq = equation.split(' <= ')
    left = eq[0].split(' + ')
    cleft = []
    for l in left:
        temp = re.search(r'(cpd\d\d\d\d\d)', l)
        item = temp.group(1)
        cleft.append(item)
    right = eq[1].split(' + ')
    left = cleft
    cright = []
    for r in right:
        temp = re.search(r'(cpd\d\d\d\d\d)', r)
        item = temp.group(1)
        cright.append(item)
    right = cright
    if '>' or '=' in direction:
        orgreacs[org][p]['substrates'] = left
        orgreacs[org][p]['products'] = right
    if '<' in direction:
        orgreacs[org][p]['substrates'] = right
        orgreacs[org][p]['products'] = left
    return


def comp2ikKegg(complist):
    tcomp = []
    for c in complist:
        if c in kegg2ik:
            tcomp.append(kegg2ik[c])
    return tcomp


def comp2ikMetacyc(complist):
    tcomp = []
    for c in complist:
        if c in metacyc2ik:
            tcomp.append(metacyc2ik[c])
    return tcomp


def comp2ikseed(complist):
    tcomp = []
    for c in complist:
        if c in cpd2ik:
            tcomp.append(cpd2ik[c])
    return tcomp


def prodScore(fblock, inchi_k, kr_list, pr_list, org):
    # calculates the production score of a compound based on deltaG values of all the reactions its in
    dgs_l = []
    reac2dg_l={}
    reactions=set()
    if org not in scores:
        scores[org]={}
    # first deal with the kegg list:
    if kr_list:
        for kreac in kr_list:
            # get the cpd id for the inchi key
            if inchi_k in ik2cpd:
                sname = ik2cpd[inchi_k]
                if  org in orgreacs and kreac in orgreacs[org]:
                    torg = orgreacs[org][kreac]
                    if sname in torg['substrates']:
                        # check if this reaction is already updated for another isomer
                        if fblock in scores[org]and 'reactions' in scores[org][fblock]:
                            if kreac in scores[org][fblock]['reactions']:
                                continue
                        dgs_l.append(torg['deltag'])
                        reac2dg_l[kreac]=torg['deltag']
                        reactions.add(kreac)
                    elif sname in torg['products']:
                        # check if this reaction is already updated for another isomer
                        if fblock in scores[org]and 'reactions' in scores[org][fblock]:
                            if kreac in scores[org][fblock]['reactions']:
                                continue
                        dgs_l.append(-1 * torg['deltag'])
                        reac2dg_l[kreac] = -1*torg['deltag']
                        reactions.add(kreac)
            elif 'kegg' in mdm[inchi_k]:
                sname = mdm[inchi_k]['kegg']
                if  org in orgreacs and kreac in orgreacs[org]:
                    torg = orgreacs[org][kreac]
                    if sname in torg['substrates']:
                        # check if this reaction is already updated for another isomer
                        if fblock in scores[org]and 'reactions' in scores[org][fblock]:
                            if kreac in scores[org][fblock]['reactions']:
                                continue
                        dgs_l.append(torg['deltag'])
                        reac2dg_l[kreac] = torg['deltag']
                        reactions.add(kreac)
                    elif sname in torg['products']:
                        # check if this reaction is already updated for another isomer
                        if fblock in scores[org]and 'reactions' in scores[org][fblock]:
                            if kreac in scores[org][fblock]['reactions']:
                                continue
                        dgs_l.append(-1 * torg['deltag'])
                        reac2dg_l[kreac] = -1*torg['deltag']
                        reactions.add(kreac)
            # tempk=sr[p]
    # plantcyc list
    if pr_list:
        for preac in pr_list:
            if 'plantcyc' in mdm[inchi_k]:
                pname = mdm[inchi_k]['plantcyc']
                if org in orgreacs and preac in orgreacs[org]:
                    tempp = orgreacs[org][preac]
                else:
                    continue
                if pname in tempp['substrates']:
                    if 'deltag' in tempp:
                        # check if this reaction is already updated for another isomer
                        if fblock in scores[org]and 'reactions' in scores[org][fblock]:
                            if preac in scores[org][fblock]['reactions']:
                                continue
                        dgs_l.append(tempp['deltag'])
                        reac2dg_l[preac] = tempp['deltag']
                        reactions.add(preac)

                elif pname in tempp['products']:
                    if 'deltag' in tempp:
                        # check if this reaction is already updated for another isomer
                        if fblock in scores[org]and 'reactions' in scores[org][fblock]:
                            if preac in scores[org][fblock]['reactions']:
                                continue
                        dgs_l.append(-1 * tempp['deltag'])
                        reac2dg_l[preac] = -1*tempp['deltag']
                        reactions.add(preac)
            else:
                # no plantcyc id how can I get the dg value
                print('no plantcyc id')

    # sum it up
    if dgs_l:
        # check I'm not repeating the exact inchi key (why??)
        if  org in visited and inchi_k in visited[org]:
            return
        if org not in scores:
            scores[org]={}
        if fblock not in scores[org]:
            scores[org][fblock]={}
        if inchi_k in ik2mass:
            scores[org][fblock]['mass']=ik2mass[inchi_k]['MW']
        if 'inchikey' not in scores[org][fblock]:
            scores[org][fblock]['inchikey']=set()
        scores[org][fblock]['inchikey'].add(inchi_k)
        if org not in visited:
            visited[org]=set()
        visited[org].add(inchi_k)
        if 'sum' not in scores[org][fblock]:
            scores[org][fblock]['sum']=0
        if 'reac_num' not in scores[org][fblock]:
            scores[org][fblock]['reac_num']=0
        if 'reactions' not in scores[org][fblock]:
            scores[org][fblock]['reactions']=set()
        if 'dgs' not in scores[org][fblock]:
            scores[org][fblock]['dgs']=[]

        scores[org][fblock]['sum'] =scores[org][fblock]['sum'] + sum(dgs_l)
        scores[org][fblock]['reac_num'] =scores[org][fblock]['reac_num']+ len(dgs_l)
        scores[org][fblock]['reactions']= scores[org][fblock]['reactions']| reactions
        scores[org][fblock]['dgs']=scores[org][fblock]['dgs']+dgs_l
        scores[org][fblock]['average'] = statistics.mean(scores[org][fblock]['dgs'])
        scores[org][fblock]['median'] = statistics.median(scores[org][fblock]['dgs'])
        '''if 'sinks' not in scores[org][fblock]:
            scores[org][fblock]['sinks']=set()'''



# ----------- main ---------

visited={}
keggCache={}
k = KEGG()
global ikscores
ikscores = {}
# load MDM
print('loading MDM')
mdm = {}
with open('MDM_021721.json') as mdmin:
    mdm = json.load(mdmin)

# load the masses for the inchi key
ik2mass=pickle.load(open('ik2fbNmass.pkl', 'rb'))


# create a mdm2kegg and mdm2metacyc convertion dictionary to validate inchi keys
global notinmdm
notinmdm = {}
global kegg2mdm
kegg2mdm = {}
global meta2mdm
meta2mdm = {}
for mik in mdm:
    if 'kegg' in mdm[mik]:
        if isinstance(mdm[mik]['kegg'], list):
            for kgg in mdm[mik]['kegg']:
                kegg2mdm[kgg] = mik
        else:
            kegg2mdm[mdm[mik]['kegg']] = mik
    if 'plantcyc' in mdm[mik]:
        meta2mdm[mdm[mik]['plantcyc']] = mik

print('loaded MDM')
print('loading seed reactions')
seedReac = {}
with open('reactions.json') as rin:
    seedReac = json.load(rin)

# change the seed dictionary to something more easy to work with
global sr
sr = {}
ccount={}
rout=open('reaction_to_check.txt','w')
for n, reac in enumerate(seedReac):
    # test: check if a reaction alias appears in more than one reaction
    if 'aliases' in reac and reac['aliases']:
        for a in reac['aliases']:
            if 'KEGG:' in a:
                place=a.split(': ')
                if place[1] in sr and reac['deltag']!= sr[place[1]]['deltag']:
                    if reac['abbreviation'] not in ccount:
                        ccount[reac['abbreviation']]=0
                    ccount[reac['abbreviation']]+=1
                    reac['abbreviation']=reac['abbreviation']+'_iso'+str(ccount[reac['abbreviation']])
                    rout.write(reac['id']+'\t'+place[1]+'\t'+ reac['definition']+ '\t'+ str(reac['deltag'])+'\n')
                    rout.write(sr[place[1]]['rxn']+'\t'+place[1]+'\t'+ sr[place[1]]['definition']+'\t'+str(sr[place[1]]['deltag'])+'\n\n')
                    print (place[1]+' -caution, multiple reactions for the same kegg rnum')
    id = reac['abbreviation']
    dg = reac['deltag']
    drc = reac['reversibility']
    definition = reac['equation']
    ndef=reac['definition']
    sr[id] = {}
    sr[id]['deltag'] = dg
    sr[id]['direction'] = drc
    sr[id]['equation'] = definition
    sr[id]['definition']=ndef
    if 'id' not in reac:
        print('no id for: '+ reac['abbreviation'])
    sr[id]['rxn']=reac['id']
rout.close()
print('loading seed compounds')
seedcomp = {}
global kegg2ik
kegg2ik = {}
global metacyc2ik
metacyc2ik = {}
global cpd2ik
cpd2ik = {}
global ik2cpd
ik2cpd = {}
global ik2kegg
ik2kegg = {}
global ik2metacyc
ik2metacyc = {}
global kegg2cpd
kegg2cpd = {}
global metacyc2cpd
metacyc2cpd = {}
global ik2name
ik2name = {}
cofactors = set()
with open('compounds.json') as sin:
    seedcomp = json.load(sin)

# turning seed compounds to a managable file
noik = 0
seediks = set()
for comp in seedcomp:
    ik = comp['inchikey'].strip()
    if ik:
        seediks.add(ik)
    if comp['is_cofactor'] == 1:
        cofactors.add(ik)
    cpd2ik[comp['id']] = ik
    ik2cpd[ik] = comp['id']
    if comp['aliases']:
        for name in comp['aliases']:
            if name.startswith('Name'):
                name1 = name.split(':')
                names = name1[1].split(';')
                ik2name[ik] = names
            if name.startswith('KEGG'):
                name1 = name.split(':')
                cs = name1[1].split(';')
                for cnum in cs:
                    cnum = cnum.strip()
                    kegg2ik[cnum] = ik
                    ik2kegg[ik] = cnum
                    kegg2cpd[cnum] = comp['id']
            if name.startswith('MetaCyc'):
                name2 = name.split(':')
                cs = name2[1].split(';')
                for cnum in cs:
                    cnum = cnum.strip()
                    metacyc2ik[cnum] = ik
                    ik2metacyc[ik] = cnum
                    metacyc2cpd[cnum] = comp['id']
    else:
        noik += 1

seediks = list(seediks)
print('loading PMN reaction data')
global pmnReac
pmnReac = {}
with open('PMNReacs.json') as pin:
    pmnReac = json.load(pin)

# step 1 - get all the reactions for an organism- in this case apple
# for each reaction check if has deltaG. If not try to get from seed.
global orgreacs
orgreacs = {}
global scores
scores = {}
totreac = set()
dgreac = set()
badkegg = 0
badpmn = 0
dgs = []
fout = open('bad_reac.txt', 'w')
pout = open('deltag.txt', 'w')
appleiks=set()
global org2syn
org2syn={}
# before starting the loop verify all inchikeys are ok
# b=verifyIk(mdm.keys(),seediks)
u = UniChem()
reac2dg={}
size_discard=set()
errs=set()
for ik in tqdm(mdm):

    '''working with first block inchi after steromers made a mess in the output'''
    fb=mdm[ik]['fblock']
    # collecting reactions for input into prodScore
    krlist = set()
    prlist = set()
    if 'orgs' in mdm[ik]:
        if ik in ik2mass and ik2mass[ik]['MW']:
            mass = ik2mass[ik]['MW']
            if not isinstance(mass, float):
                if isinstance(mass,str) and re.match(r"\d+\.\d+[a-zA-Z]\d+",mass):
                    mass = mass[:-2]
                mass=float(mass)
            if mass < 50:
                size_discard.add(ik)
                continue
        else:
            errs.add(ik)
            continue
        for org in mdm[ik]['orgs']:
            mdmorg=None
            pmnorg=None
            # insert an organism verification mechanism
            if 'subsp.' in org:
                temp = org.split('subsp.')
            elif 'var.' in org:
                temp=org.split('var.')
            else:
                temp = org.split('subsp.')
            temp[0] = temp[0].strip()
            if temp[0] in pmnReac:
                org2syn[org] = temp[0]
                mdmorg = org
                pmnorg = temp[0]
            #appleiks.add(ik)

            for r in mdm[ik]['orgs'][org]:
                if not re.match(r'\d+\.\d+\.\d+\.\d+$',r):
                    # plantcyc reaction
                    if 'Zea mays' in org:
                        pmnorg='Zea mays mays'
                        mdmorg='Zea mays'
                    if 'Citrullus lanatus' in org:
                        pmnorg='Citrullus lanatus vulgaris'
                        mdmorg='Citrullus lanatus'
                    if 'Oryza sativa japonica'in org:
                        pmnorg='Oryza sativa Japonica Group'
                        mdmorg='Oryza sativa japonica'
                    if 'Camellia sinensis' in org:
                        pmnorg='Camellia sinensis sinensis'
                        mdmorg='Camellia sinensis'
                    if 'Beta vulgaris' in org:
                        pmnorg = 'Beta vulgaris vulgaris'
                        mdmorg = 'Beta vulgaris'
                    if 'Olea europaea var. sylvestris' in org:
                        pmnorg='Olea europaea'
                        mdmorg='Olea europaea var. sylvestris'
                    if 'Hordeum vulgare subsp. Vulgare' in org:
                        pmnorg='Hordeum vulgare'
                        mdmorg='Hordeum vulgare subsp. Vulgare'
                    if pmnorg:
                        org=pmnorg
                    if r in pmnReac[org]:
                        if mdmorg not in orgreacs:
                            orgreacs[mdmorg]={}
                        orgreacs[mdmorg][r] = {}
                        prlist.add(r)
                        totreac.add(r)
                        if 'deltag' in pmnReac[org][r]:
                            orgreacs[mdmorg][r]['deltag'] = pmnReac[org][r]['deltag']
                            reac2dg[r]=pmnReac[org][r]['deltag']
                            dgreac.add(r)
                            '''pout.write(
                                r + '\t' + str(pmnReac[org][r]['deltag']) + '\t' + ','.join(pa.keys()) + '\n')'''
                            dgs.append(pmnReac[org][r]['deltag'])
                        orgreacs[mdmorg][r]['substrates'] = pmnReac[org][r]['substrates']
                        orgreacs[mdmorg][r]['products'] = pmnReac[org][r]['products']
                    else:
                        #print('what now?')
                        badpmn += 1
                else:
                    if mdmorg:
                        org=mdmorg
                    if 'Zea mays mays' in org:
                        org = 'Zea mays'
                    if 'Citrullus lanatus vulgaris' in org:
                        org='Citrullus lanatus'
                    if 'Oryza sativa Japonica Group' in org:
                        org='Oryza sativa japonica'
                    if 'Beta vulgaris vulgaris' in org:
                        org='Beta vulgaris'
                    if 'Camellia sinensis sinensis' in org:
                        org='Camellia sinensis'
                    if 'Olea europaea' in org:
                        org='Olea europaea var. sylvestris'
                    if 'Hordeum vulgare' in org:
                        org='Hordeum vulgare subsp. Vulgare'

                    pa=mdm[ik]['orgs'][org][r]
                    for p in pa:
                        if p in sr:
                            totreac.add(p)
                            if org not in orgreacs:
                                orgreacs[org] = {}
                            # a kegg situation, r is in fact an enzyme
                            orgreacs[org][p] = {}
                            krlist.add(p)
                            if p in keggCache:
                                kreac=keggCache[p]
                            else:
                                kreac = k.parse(k.get(p))
                                keggCache[p]=kreac
                            if 'deltag' in sr[p]:
                                if sr[p]['deltag'] == 10000000.0:
                                    fout.write(p + '\t' + sr[p]['equation'] + '\n')
                                    badkegg += 1
                                    del orgreacs[org][p]
                                    krlist.remove(p)
                                    continue
                                orgreacs[org][p]['deltag'] = sr[p]['deltag']
                                reac2dg[p] = sr[p]['deltag']
                                dgreac.add(p)
                                pout.write(p + '\t' + str(sr[p]['deltag']) + '\t' + ','.join(pa[p]) + '\n')
                                dgs.append(sr[p]['deltag'])
                            if isinstance(kreac, dict) and 'EQUATION' in kreac:
                                getSubsProdsKEGG(kreac['EQUATION'], sr[p]['direction'],org)
                            else:
                                getSubsProdsSEED(sr[p]['equation'], sr[p]['direction'],org)
                        else:
                            #print(p + ' bummer')
                            '''if isinstance(kreac,dict):
                                fout.write(p + '\t' + kreac['EQUATION'] + '\n')'''
                            badkegg += 1

            if mdmorg:
                org=mdmorg
            if 'Zea mays mays' in org:
                org = 'Zea mays'
            if 'Citrullus lanatus vulgaris' in org:
                org='Citrullus lanatus'
            if 'Oryza sativa Japonica Group' in org:
                org='Oryza sativa japonica'
            if 'Beta vulgaris vulgaris' in org:
                org='Beta vulgaris'
            if 'Camellia sinensis sinensis' in org:
                org='Camellia sinensis'
            if 'Olea europaea' in org:
                org='Olea europaea var. sylvestris'
            if 'Hordeum vulgare' in org:
                org='Hordeum vulgare subsp. Vulgare'
            if 'Brassica rapa' in org:
                org='Brassica rapa subsp. Pekinensis'
            # finished collecting both types of reactions
            if not krlist and not prlist:
                continue
            else:
                prodScore(fb, ik, krlist, prlist,org)
# save orgreacs to file

pickle.dump( orgreacs, open( "orgReacKinetics.pkl", "wb" ) )
# save scores to file
# change organism names to lower for scores.
for org in scores:
    norg=org.lower()
    scores[norg]=scores.pop(org)
pickle.dump( scores, open( "kineticScores.pkl", "wb" ) )

print('total reacs ' + str(len(orgreacs)) + ' reactions')
print('total compounds '+ str(len(mdm.keys())))
print('total scored compounds'+ str(len(scores.keys())))

venn2([totreac, dgreac], set_labels=('Total', 'Has DG'), set_colors=('steelblue', 'skyblue'), alpha=0.7)

plt.show()

# deltag distribution
df = pd.DataFrame({'Delta G': dgs})
plt.hist(dgs, bins=50, facecolor='blue', alpha=0.5)
plt.show()
sns.distplot(df, bins=50, hist=True, kde=True, color='steelblue',
             hist_kws={'edgecolor': 'darkblue'},
             kde_kws={'linewidth': 1})
plt.show()
'''# scores distribution
scdf=pd.DataFrame.from_dict(scores, orient='index')
sns.distplot(df, bins=50, hist=True, kde=True, color='steelblue',
             hist_kws={'edgecolor': 'darkblue'},
             kde_kws={'linewidth': 1})
sns.distplot(scdf['sum'],bins=50, hist=True, kde=True, color = 'mediumseagreen',
             hist_kws={'edgecolor':'seagreen'},
             kde_kws={'linewidth': 1})
plt.show()
sns.distplot(scdf['reac_num'],bins=50, hist=True, kde=True, color = 'gold',
             hist_kws={'edgecolor':'orange'},
             kde_kws={'linewidth': 1})
plt.show()
sns.distplot(scdf['sinknum'],bins=50, hist=True, kde=True, color = 'crimson',
             hist_kws={'edgecolor':'firebrick'},
             kde_kws={'linewidth': 1})
plt.show()'''
# output the scores
sout = open('dg_scores.txt', 'w')
sout.write('Inchi key\tName\tSum dg\tReaction number\tAverage dg\tMedian dg\tpotential sinks\tReactions and dg\n')
for g,s in enumerate(scores):
    print('score: '+ str(g))
    if 'name' in mdm[s]:
        if 'sinks' in scores[s]:
            sout.write(s + '\t' + mdm[s]['name'] + '\t' + str(scores[s]['sum']) + '\t' + str(scores[s]['reac_num']) + '\t'
                       + str(scores[s]['average']) + '\t' + str(scores[s]['median']) +'\t'+
                       ','.join(list(scores[s]['sinks'])) +'\t'+str(scores[s]['sinknum'])+
                       '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
        else:
            sout.write(
                s + '\t' + mdm[s]['name'] + '\t' + str(scores[s]['sum']) + '\t' + str(scores[s]['reac_num']) + '\t'
                + str(scores[s]['average']) + '\t' + str(scores[s]['median']) + '\tno sinks\tNA'+
                '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
    else:
        while True:
            try:
                pcomp=pcp.get_compounds(s,'inchikey')
                if pcomp and hasattr(pcomp[0],'synonyms') and pcomp[0].synonyms:
                    if 'sinks' in scores[s]:
                        sout.write(
                            s + '\t'+pcomp[0].synonyms[0]+'\t' + str(scores[s]['sum']) + '\t' + str(scores[s]['reac_num']) + '\t'
                            + str(scores[s]['average']) + '\t' + str(scores[s]['median']) +'\t'+
                            ','.join(list(scores[s]['sinks']))+'\t'+str(scores[s]['sinknum'])+
                            '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
                        break
                    else:
                        sout.write(
                            s + '\t' + pcomp[0].synonyms[0] + '\t' + str(scores[s]['sum']) + '\t' + str(
                                scores[s]['reac_num']) + '\t'
                            + str(scores[s]['average']) + '\t' + str(scores[s]['median']) + '\tno sinks\tNA'+
                            '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
                        break
                elif pcomp and hasattr(pcomp[0],'name') and pcomp[0].name:
                    if 'sinks' in scores[s]:
                        sout.write(
                            s + '\t'+pcomp[0].name+'\t' + str(scores[s]['sum']) + '\t' + str(scores[s]['reac_num']) + '\t'
                            + str(scores[s]['average']) + '\t' + str(scores[s]['median']) +'\t'+
                            ','.join(list(scores[s]['sinks']))+ '\t'+str(scores[s]['sinknum'])+
                            '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
                        break
                    else:
                        sout.write(
                            s + '\t' + pcomp[0].name + '\t' + str(scores[s]['sum']) + '\t' + str(
                                scores[s]['reac_num']) + '\t'
                            + str(scores[s]['average']) + '\t' + str(scores[s]['median']) + '\tno sinks\tNA'+
                            '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
                        break
                else:
                    if 'sinks' in scores[s]:
                        sout.write(
                            s + '\tInsert name here\t' + str(scores[s]['sum']) + '\t' + str(
                                scores[s]['reac_num']) + '\t'
                            + str(scores[s]['average']) + '\t' + str(scores[s]['median']) +'\t'+
                            ','.join(list(scores[s]['sinks']))+ '\t'+str(scores[s]['sinknum'])+
                            '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
                        break
                    else:
                        sout.write(
                            s + '\tInsert name here\t' + str(scores[s]['sum']) + '\t' + str(
                                scores[s]['reac_num']) + '\t'
                            + str(scores[s]['average']) + '\t' + str(scores[s]['median']) + '\tno sinks\tNA'+
                            '\t'+ json.dumps(scores[s]['reac2dg'])+'\n')
                        break
            except Exception as ex:
                # If an error is raised it will print it while the code runs.
                exc_type, exc_obj, exc_tb = sys.exc_info()
                print(exc_type, exc_tb.tb_lineno)
                print('Retrying in 0.5 seconds. Exception in connecting')
                sleep(0.5)


fout.close()
pout.close()
sout.close()

'''# get experimental inchi keys and get the subset of those with scores
# read the experimental platform file
appleexpkeys = set()
with open('exp_platform.txt') as ein:
    for line in ein:
        if line.startswith('Name'):
            continue
        line = line.split('\t')
        expinchikey = line[1].strip()
        if '1' in line[3] or '1' in line[8] or '1' in line[13] or '1' in line[18] or '1' in line[23]:
            appleexpkeys.add(expinchikey)

# compare scores to experimental
skeys=list(scores.keys())
venn2([set(skeys), appleexpkeys], set_labels=('Have score', 'Experimental'), set_colors=('darkseagreen', 'seagreen'), alpha=0.7)
plt.show()
tout=open('exp_dg_score_apple.txt','w')
for ap in appleexpkeys:
    if ap in scores:
        if 'sinks' in scores[ap] and 'name' in mdm[ap]:
            tout.write(
                ap +  '\t' + mdm[ap]['name']+ '\t'+str(scores[ap]['sum']) + '\t' + str(scores[ap]['reac_num']) + '\t'
                + str(scores[ap]['average']) + '\t' + str(scores[ap]['median']) + '\t' +
                ','.join(list(scores[ap]['sinks'])) + '\t' + str(scores[ap]['sinknum']) + '\n')
        else:
            tout.write(
                ap + '\t No name\t' + str(scores[ap]['sum']) + '\t' + str(scores[ap]['reac_num']) + '\t'
                + str(scores[ap]['average']) + '\t' + str(scores[ap]['median']) + '\tno sinks\tNA\n')
tout.close()'''

with open('reac2dg.json','w') as nin:
    nin.write(json.dumps(reac2dg))
sys.exit()
