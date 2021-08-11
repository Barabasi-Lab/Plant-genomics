''' analyze the metabolism in mdm according to plant familis and visualize.
1. number of plants per family - what's the average
2. unique metabolites per family and unique to family metabolites
3. primary/secondary metabolism from a family pov
Written by Shany Ofaim, the Barabasi lab CCNR Northeastern University 2021'''


# imports
import json
import statistics
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np
from ete3 import NCBITaxa
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
from ete3 import PhyloTree
import sys
from ete3 import Nexml

# functions
def mylayout(node):
    # If node is a leaf, add the nodes name and a its scientific
    node.name=node.sci_name
    # name
    if node.is_leaf():
        N = AttrFace("sci_name", fsize=35)
        faces.add_face_to_node(N, node, 0, position="aligned")

# main

# classify organisms to families
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

orgList=set()
org2family={}
family2orgs={}
with open('C:\Barabasi\Foodome\Plant_genomes\KEGG_foodlist.txt') as fin: #your organism list goes here
    for line in fin:
        line=line.split('\t')
        name=line[0].strip()
        family=line[1].strip()
        orgList.add(name)
        org2family[name]=family
        if family not in family2orgs:
            family2orgs[family]=set()
        family2orgs[family].add(name)

# load mdm
print('loading MDM')
mdm={}
with open('MDM.json') as mdmin:
    mdm=json.load(mdmin)

# load primary/secondary definitions
# add primary and secondary info
cat2class={}
with open('Primary_Secondary.txt') as psin:
    for line in psin:
        line=line.split('\t')
        pscat=line[0].strip()
        psclass=line[1].strip()
        cat2class[pscat]=psclass
# read pathway codes
pathways={}

with open('pathwayCodes.txt') as pwin:
    for line in pwin:
        if line is '':
            continue
        line=line.split('\t')
        pcode=line[0].strip()
        pname=line[1].strip()
        pcat=line[2].strip()
        if pname not in pathways:
            pathways[pname]={}
        pathways[pname]['cat']=pcat
        pathways[pname]['code']=pcode
        if pcat:
            pathways[pname]['class']=cat2class[pcat]

# start a metabolite centric dictionary
compounds={}
orgs=set()
orgnames=set()
notfound=set()
taxid2org={}
for ik in mdm:
    # only looking at compounds with organisms
    if 'orgs' in mdm[ik]:
        comp2families=set()
        for org in mdm[ik]['orgs']:
            orgnames.add(org)
            name2taxid = ncbi.get_name_translator([org])
            if name2taxid:
                orgs.add(name2taxid[org][0])
                taxid2org[name2taxid[org][0]]=org
            else:
                notfound.add(org)
            if ik not in compounds:
                compounds[ik]={}
            # get pathways
            '''for r in mdm[ik]['orgs'][org]:
                if r in mdm[ik]['orgs'][org]:
                    pa = mdm[ik]['orgs'][org][r]
                    for p in pa:
                        if not isinstance(pa[p], dict):
                            # plantcyc
                            for path in pa:
                                print (pa)
                        else:
                            # a kegg situation, r is in fact an enzyme

                            for path in pa[p]:
                                print(pa[p])'''

fout=open('taxid2org.txt','w')
for tax in taxid2org:
    fout.write(str(tax)+'\t'+taxid2org[tax]+'\n')
fout.close()
tree = ncbi.get_topology(list(orgs))
tree1 = PhyloTree('(((((((4565:1,112509:1)1:1,(4528:1,4538:1,4537:1,65489:1)1:1,38705:1)1:1,((4558:1,4577:1)1:1,4555:1)1:1)1:1,4615:1)1:1,(42345:1,51953:1)1:1,(4641:1,52838:1)1:1)1:1,55577:1,4686:1)1:1,((((((((3716:1)1:1,(51351:1)1:1,3708:1)1:1,3726:1)1:1,90675:1)1:1,3649:1)1:1,((2711:1,85681:1)1:1,171929:1)1:1,(3641:1,66656:1)1:1)1:1,((((((23211:1,225117:1)1:1,3750:1)1:1,(3760:1,102107:1,42229:1)1:1)1:1,57918:1)1:1,3483:1,326968:1)1:1,(3983:1,3988:1)1:1,(((3656:1,3659:1)1:1,3654:1)1:1,(3661:1,3662:1,3664:1)1:1,3673:1)1:1,((((3914:1,3917:1,157791:1)1:1,3885:1,3821:1,3847:1)1:1,3827:1)1:1,3871:1)1:1,51240:1)1:1,29760:1)1:1,(((4232:1,4236:1,59895:1)1:1,4039:1)1:1,((4182:1,226208:1,158386:1)1:1,((4113:1,4111:1,(28526:1,4081:1)1:1)1:1,4072:1)1:1)1:1,4442:1)1:1,((63459:1,3562:1)1:1,161934:1)1:1)1:1,4432:1);', sp_naming_function=lambda name: name)
nout=open('node_names.txt','w')
distMtx= {}
distclust={}
for node1 in tree.traverse("levelorder"):
  # Do some analysis on node
  nout.write(node1.sci_name+'\n')
  node1.name=node1.sci_name
  print (node1.name)
  if 'Dioscorea cayenensis subsp. rotundata' in node1.name:
      node1.name='Dioscorea rotundata'
  if node1.name in orgnames:
      if node1.name not in distMtx:
          distMtx[node1.name]={}
      for node2 in tree.traverse("postorder"):
          if 'Dioscorea cayenensis subsp. rotundata' in node2.sci_name:
              node2.sci_name='Dioscorea rotundata'
          if node2.sci_name in orgnames:

              dist=tree.get_distance(node1,node2)
              if dist not in distclust:
                  distclust[dist]=set()
              distclust[dist].add(node1.name)
              distclust[dist].add(node2.sci_name)
              distMtx[node1.name][node2.sci_name] = dist
              #print('distance '+ node1.name+ ' '+ node2.sci_name+' is: '+str(dist))
nout.close()
df=pd.DataFrame(distMtx)

# get common ancestors for the nodes on the tree for the family classification
#node2labels = tree.get_cached_content(store_attr=["name",'sci_name','rank'])
fnode=tree.search_nodes(rank='clade')
family2members={}
ancestores=[]
for fn in fnode:
    name=fn.sci_name
    ancestores.append(name)
    if name not in family2members:
        family2members[name]={}
        family2members[name]['snodes'] = []
        family2members[name]['nodes']=[]
    # get the members
        for f in fn:
            nname=f.name
            sname=f.sci_name
            family2members[name]['snodes'].append(sname)
            family2members[name]['nodes'].append(nname)


for n in tree.traverse():
   nstyle = NodeStyle()
   nstyle["hz_line_width"] = 3
   nstyle["vt_line_width"] = 3
   nstyle["vt_line_color"] = "#909497"
   nstyle["hz_line_color"] = "#909497"
   n.set_style(nstyle)
circular_style = TreeStyle()
#circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 100
circular_style.rotation = 90
circular_style.layout_fn=mylayout
circular_style.show_leaf_name=False
circular_style.guiding_lines_color="#909497"
circular_style.draw_guiding_lines=True
circular_style.layout_fn=mylayout
circular_style.show_leaf_name=False
tree.render("mytree.png", w=5000, units="mm", tree_style=circular_style)
tree.show()
#tree.set_species_naming_function()
tree.write(format=1, outfile="new_tree.nw")
print (tree.get_ascii(attributes=['sci_name', 'taxid','rank']))



