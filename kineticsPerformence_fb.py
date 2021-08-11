''' Reviewing the kinetics performance in different rankings windows then comparing to the same process for
the unranked genome vs. experment
Written by Shany Ofaim, the Barabasi lab CCNR northeastern university 2021'''


# imports
import pickle
from tqdm import tqdm
import seaborn as sns
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from bioservices import *
from sklearn.metrics import confusion_matrix
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import plot_precision_recall_curve
from numpy import sqrt
from numpy import argmax
from sklearn.metrics import f1_score
from sklearn.metrics import auc


# ------------ functions -------------
def convert2fb(iklist):
    ''' converts a full key list to first block list'''
    fblist=set()
    for ik in iklist:
        ikl=ik.split('-')
        fblist.add(ikl[0])
        fb2ik[ikl[0]]=ik
    return fblist





# ----------- main ------------------
fb2ik={}
''' This analysis will be performed for each genome of the 13 that we have experiments for '''
''' edit 06/14/21 - add the complete genome without kinetic scores and calculate the perfomrance of guessing the experiments'''
# load the file with the kinetics scores - contains both positives and negatives
orgdgscores=pickle.load(open('kineticScores.pkl','rb'))
for org in orgdgscores:
    norg=org.lower()
    orgdgscores[norg]=orgdgscores.pop(org)
# load the file with all the keys from different sources including experimental - full keys
orgkeys=pickle.load(open('orgkeys_new.pkl','rb'))

# get the list of plants that have experimental evidence
plantList=set()
org2total={}
mtdf=pd.DataFrame()
# get all known fbs
totfbs=set()
with open('all_fb.txt') as fin:
    for line in fin:
        line=line.strip()
        totfbs.add(line)

for org in tqdm(orgkeys):
    if 'experiment' in orgkeys[org]:
        plantList.add(org)
        org2total[org]=len(orgkeys[org]['total'])
        # get the genomics fraction:
        # genome set
        if 'kegg' not in orgkeys[org]:
            gset = orgkeys[org]['pmn']
        if 'pmn' not in orgkeys[org]:
            gset = orgkeys[org]['kegg']
        if 'kegg' in orgkeys[org] and 'pmn' in orgkeys[org]:
            gset = orgkeys[org]['kegg'] | orgkeys[org]['pmn']
        # manually curated inchi key discrepencies with NDM fix as a result of keys describing the same molecule but coming from different sources (unichem, PMN)
        cur = ['LFTYTUAZOPRMMI-MPIASZHXSA-N', 'SRBFZHDQGSBBOR-SOOFDHNKSA-N', 'MEFKEPWMEQBLKI-XCPQSEKJSA-O',
               'BAWFJGJZGIEFAR-NNYOXOHSSA-N', 'CBIDVWSRUUODHL-QTSLKERKSA-N', 'VOXXWSYKYCBWHO-UHFFFAOYSA-N',
               'JPIJQSOTBSSVTP-STHAYSLISA-M', 'OSJPPGNTCRNQQC-UHFFFAOYSA-N', 'NKDFYOWSKOHCCO-YPVLXUMRSA-N',
               'JZNWSCPGTDBMEW-UHFFFAOYSA-N', 'SBLKVIQSIHEQOF-UPHRSURJSA-N', 'RBNPOMFGQQGHHO-UHFFFAOYSA-N',
               'JYJIGFIDKWBXDU-MNNPPOADSA-N', 'FYSSBMZUBSBFJL-UHFFFAOYSA-N', 'ACWASDPGAVYCNI-JKUQZMGJSA-N',
               'BRGMHAYQAZFZDJ-PVFLNQBWSA-N', 'LLCSXHMJULHSJN-UHFFFAOYSA-N', 'ODBLHEXUDAPZAU-ZAFYKAAXSA-N',
               'HEBKCHPVOIAQTA-NGQZWQHPSA-N', 'RHGKLRLOHDJJDR-UHFFFAOYSA-N', 'HXEACLLIILLPRG-UHFFFAOYSA-N',
               'ZZLHPCSGGOGHFW-ZMQIUWNVSA-N', 'VCWMRQDBPZKXKG-SPBUTQSFSA-N', 'JFCQEDHGNNZCLN-UHFFFAOYSA-N',
               'BJEPYKJPYRNKOW-UHFFFAOYSA-N', 'LKDRXBCSQODPBY-VRPWFDPXSA-N', 'HSCJRCZFDFQWRP-ABVWGUQPSA-N',
               'AAWZDTNXLSGCEK-LNVDRNJUSA-N', 'ILGMGHZPXRDCCS-UHFFFAOYSA-N', 'NGSWKAQJJWESNS-UHFFFAOYSA-N',
               'FEWJPZIEWOKRBE-UHFFFAOYSA-N', 'PADQINQHPQKXNL-UHFFFAOYSA-N', 'FDHFJXKRMIVNCQ-OSIZZBRKSA-N',
               'QQHJDPROMQRDLA-UHFFFAOYSA-N', 'QAIPRVGONGVQAS-RQOWECAXSA-N', 'YDBYJHTYSHBBAU-UHFFFAOYSA-O',
               'POJWUDADGALRAB-UHFFFAOYSA-N', 'CUZKLRTTYZOCSD-UHFFFAOYSA-N', 'BJRNKVDFDLYUGJ-UHFFFAOYSA-N']
        gset = gset | set(cur)
        # convert to first block
        gset = convert2fb(gset)
        orgkeys[org]['experiment_fb']=set()
        for k in orgkeys[org]['experiment']:
            kl=k.split('-')
            kfb=kl[0]
            orgkeys[org]['experiment_fb'].add(kfb)

        for fik in totfbs:
            if fik in gset:
                genome=1
            else:
                genome=0
            if fik in orgkeys[org]['experiment_fb']:
                experi=1
            else:
                experi=0
            # update the dataframe
            values_to_add = {'plant': org, 'genome':genome, 'experiment':experi}
            row_to_add = pd.Series(values_to_add, name=fik)
            mtdf = mtdf.append(row_to_add)




mtresdf=pd.DataFrame()
# calculate the performance of just the genome against experiment
for plant in plantList:
    plantdf = mtdf[mtdf['plant'] == plant]
    # calculate roc curve
    fpr, tpr, thresholds = roc_curve(plantdf['genome'], plantdf['experiment'])
    aucs = roc_auc_score(plantdf['experiment'], plantdf['genome'])
    # precision recall
    average_precision = average_precision_score(plantdf['experiment'], plantdf['genome'])
    precision, recall, pthresholds = precision_recall_curve(plantdf['experiment'], plantdf['genome'])
    lr_f1, lr_auc = f1_score(plantdf['experiment'], plantdf['genome']), auc(recall, precision)

    # update the dataframe
    values_to_add = {'plant': plant, 'roc_auc': aucs,
                     'pr_auc': lr_auc,
                     'f1': lr_f1, 'average_precision': average_precision}
    row_to_add = pd.Series(values_to_add, name=plant)
    mtresdf = mtresdf.append(row_to_add)

print(mtresdf.mean())
print(mtresdf.std())
gendf=pd.DataFrame()
'''for org in tqdm(orgkeys):
    # create the dataframe for performance calculations
    if org in orgdgscores:
        for ik in orgdgscores[org]:
            score=orgdgscores[org][ik]['sum']
            rnum=orgdgscores[org][ik]['reac_num']
            avg=orgdgscores[org][ik]['average']
            mass=orgdgscores[org][ik]['mass']
            fullkeys=orgdgscores[org][ik]['inchikey']


            # is it in the experimental set?
            if 'experiment' in orgkeys[org]:
                if ik in orgkeys[org]['experiment_fb']:
                    exp=1
                else:
                    exp=0
                if score>0:
                    sumcls=1
                else:
                    sumcls=0
                # update the dataframes
                # general df
                values_to_add = {'inchi key':ik,'sum':score,'no_reactions':rnum,'average':avg,
                                 'mass': mass,'experiment':exp,'plant':org, 'sum_class':sumcls,'inchi_keys':fullkeys}
                row_to_add = pd.Series(values_to_add, name=org)
                gendf = gendf.append(row_to_add)

gendf.to_csv('kinetics_all_Plants_df_fb.csv')'''
gendf=pd.read_csv('kinetics_all_Plants_df_fb.csv')
# plot the scores vs. number of reactions
'''fig, ax = plt.subplots()
ax2 = ax.twinx()
sns.regplot(x="no_reactions", y="sum", data=gendf, order=2, ax=ax, color='black')
sns.regplot(x="no_reactions", y="average", data=gendf, order=2, ax=ax2, color='dodgerblue')
ax.set_yscale('log')
ax2.set_yscale('log')

ax2.legend(handles=[a.lines[0] for a in [ax,ax2]],
           labels=["sum", "average"])
plt.tight_layout()
plt.show()'''
# plot the scre distribution to establish a cutoff
'''plt.figure(figsize=(15, 10))
sns.set(font_scale=2.5)
g1=sns.distplot(gendf['average'], color='forestgreen',kde_kws=dict(linewidth=5), label='scores')

plt.tight_layout()
plt.show()'''
# scatter only reaction number and sum
'''plt.figure(figsize=(15, 10))
#sns.scatterplot(x='no_reactions',y='sum', data=gendf, color='grey', label='sum')
sns.scatterplot(x='no_reactions',y='average', data=gendf, color='red', label='average')
plt.legend()
plt.tight_layout()
plt.show()'''

# start per genome analysis
metdf=pd.DataFrame()
for plant in tqdm(sorted(plantList)):
    if 'pyrus communis' in plant:
        continue
    print('\n'+plant)
    plantdf = gendf[gendf['plant'] == plant]
    eor=len(plantdf)
    # look for the best metrics
    for x in range(10,eor,50):
        plantdf = gendf[gendf['plant'] == plant]
        plantdf = plantdf.sort_values(by=['sum'], ascending=False)
        plantdf=plantdf.head(n=x)
        # calculate roc curve
        fpr, tpr, thresholds = roc_curve(plantdf['experiment'],plantdf['sum'])

        # calculate the geometric mean for each threshold
        gmeans = sqrt(tpr * (1 - fpr))
        # locate the index of the largest g-mean
        ix = argmax(gmeans)
        #print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))

        # calculate AUC
        try:
            aucs = roc_auc_score(plantdf['experiment'],plantdf['sum'])
        except:
            pass



        '''# plot the roc curve
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(tpr, fpr)
        ax.set_xlabel('True positive rate')
        ax.set_ylabel('False positive rate')
        ax.legend(loc='center left')
        plt.show()'''

        #precision recall
        average_precision = average_precision_score(plantdf['experiment'],plantdf['sum_class'])
        precision, recall, pthresholds = precision_recall_curve(plantdf['experiment'],plantdf['sum_class'])
        lr_f1, lr_auc = f1_score(plantdf['experiment'], plantdf['sum_class']), auc(recall, precision)
        '''fig1, ax = plt.subplots(figsize=(6, 6))
        ax.plot(recall, precision)
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.legend(loc='center left')
        plt.show()
    
        # F-1 score
        fscore=f1_score(plantdf['experiment'],plantdf['sum_class'])'''

        # update the dataframe
        values_to_add = {'plant':plant,'df_fraction':x,'best_threshold':thresholds[ix],'roc_auc':aucs,'pr_auc':lr_auc,
                         'f1':lr_f1,'average_precision':average_precision,'norm_fraction':round((x/eor),2),
                         'total_scored':eor,'total_plant':org2total[plant]}
        row_to_add = pd.Series(values_to_add, name=plant)
        metdf = metdf.append(row_to_add)

# scatter x vs auc and x vs average precision anf f1
metdf.to_csv('performance_results_fb.csv')
sns.set(font_scale=2)
sns.color_palette("pastel")
plt.figure(figsize=(10, 10))
sns.scatterplot(x='norm_fraction',y='roc_auc', data=metdf, hue='plant',palette=sns.color_palette('Paired', n_colors=12),s=60)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., fontsize='small')
plt.ylim(0,1)
plt.axhline(y=0.5, color='black', linestyle='--')
plt.tight_layout()
plt.savefig("roc_auc.svg")
plt.show()

plt.figure(figsize=(10, 10))
sns.scatterplot(x='norm_fraction',y='pr_auc', data=metdf, hue='plant',palette=sns.color_palette('Paired', n_colors=12),s=60)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., fontsize='small')
plt.ylim(0,1)
plt.axhline(y=0.5, color='black', linestyle='--')
plt.tight_layout()
plt.savefig("pr_auc.svg")
plt.show()

plt.figure(figsize=(10, 10))
sns.scatterplot(x='df_fraction',y='average_precision', data=metdf, hue='plant',palette=sns.color_palette('Paired', n_colors=12))
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., fontsize='medium')
plt.tight_layout()
plt.savefig("avg_precision.svg")
plt.show()

plt.figure(figsize=(10, 10))
sns.scatterplot(x='df_fraction',y='f1', data=metdf, hue='plant',palette=sns.color_palette('Paired', n_colors=12))
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., fontsize='medium')
plt.tight_layout()
plt.savefig("f1.svg")
plt.show()

plt.figure(figsize=(10, 10))
sns.scatterplot(x='df_fraction',y='best_threshold', data=metdf, hue='plant',palette=sns.color_palette('Paired', n_colors=12))
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., fontsize='medium')
plt.tight_layout()
plt.savefig("best_th.svg")
plt.show()



# optimal fractions
# second pass classification

othresholds=pd.read_csv('optimal_thresholds_fb.csv')

tdict=othresholds.set_index('plant').to_dict('index')

optdf=pd.DataFrame()

optorg2keys={}
for org in tqdm(orgkeys):
    if 'pyrus communis' in org:
        continue
    # create the dataframe for performance calculations
    if org in orgdgscores:

        expset=set()
        predset=set()
        genomeset=set()
        for ik in orgdgscores[org]:
            genomeset.add(ik)
            score=orgdgscores[org][ik]['sum']
            rnum=orgdgscores[org][ik]['reac_num']
            avg=orgdgscores[org][ik]['average']
            mass=orgdgscores[org][ik]['mass']
            fullkeys = orgdgscores[org][ik]['inchikey']

            # is it in the experimental set?
            if 'experiment' in orgkeys[org]:
                if ik in orgkeys[org]['experiment_fb']:
                    exp=1
                    expset.add(ik)
                else:
                    exp=0
                if score>tdict[org]['optimal_th']:
                    sumcls=1
                    predset.add(ik)
                else:
                    sumcls=0
                # update the dataframes
                # general df
                values_to_add = {'inchi_key':ik,'sum':score,'no_reactions':rnum,'average':avg,
                                 'mass': mass,'experiment':exp,'plant':org, 'sum_class':sumcls,'inchi_keys':fullkeys}
                row_to_add = pd.Series(values_to_add, name=org)
                optdf = optdf.append(row_to_add)




pal = sns.color_palette('Paired',12)
pal=pal.as_hex()
# start optimal per genome analysis
ometdf = pd.DataFrame()
ccounter=0
plt.figure(figsize=(30, 30))

for plant in tqdm(sorted(plantList)):
    if 'pyrus communis' in plant:
        continue
    print('\n' + plant)
    plantdf = optdf[optdf['plant'] == plant]
    plantdf = plantdf.sort_values(by=['sum'], ascending=False)
    plantdf = plantdf.head(n=tdict[plant]['fraction'])
    fname=plant+'_top fraction_avg.csv'
    plantdf.to_csv(fname)

    # turn the dataframe to a dictionary to get the set of full inchi keys per optimal fraction
    plantdict=pd.DataFrame.to_dict(plantdf,orient='records')
    fullkeys=set()
    for item in plantdict:
        fullkeys=fullkeys|item['inchi_keys']

    optorg2keys[plant]=fullkeys
    # calculate roc curve
    fpr, tpr, thresholds = roc_curve(plantdf['experiment'], plantdf['sum'])

    # calculate the geometric mean for each threshold
    gmeans = sqrt(tpr * (1 - fpr))
    # locate the index of the largest g-mean
    ix = argmax(gmeans)
    # print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))

    # calculate AUC

    aucs = roc_auc_score(plantdf['experiment'], plantdf['sum'])

    # plot the roc curve
    sns.set(font_scale=3)
    plt.plot(fpr, tpr, label=plant+', auc='+str(round(aucs,2)),
             color=pal[ccounter],linewidth=4)
    plt.xlim([0, 1])
    plt.ylim([0, 1])

    ccounter+=1

    '''# precision recall
    average_precision = average_precision_score(plantdf['experiment'], plantdf['sum_class'])
    precision, recall, pthresholds = precision_recall_curve(plantdf['experiment'], plantdf['sum_class'])
    lr_f1, lr_auc = f1_score(plantdf['experiment'], plantdf['sum_class']), auc(recall, precision)
    ''''''fig1, ax = plt.subplots(figsize=(6, 6))
    ax.plot(recall, precision)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.legend(loc='center left')
    plt.show()

    # F-1 score
    fscore=f1_score(plantdf['experiment'],plantdf['sum_class'])''''''

    # update the dataframe
    values_to_add = {'plant': plant, 'df_fraction': tdict[plant]['fraction'], 'best_threshold': thresholds[ix], 'roc_auc': aucs,
                     'pr_auc': lr_auc,
                     'f1': lr_f1, 'average_precision': average_precision}
    row_to_add = pd.Series(values_to_add, name=plant)
    metdf = metdf.append(row_to_add)'''
plt.axline([0, 0], [1, 1], linewidth=6,color='grey',linestyle='--' )
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0., fontsize=28)
plt.savefig('opt_ROCs.svg')
plt.show()

pickle.dump( optorg2keys, open( "optimal_orgs.pkl", "wb" ) )

print('Done')
