import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from mol2vec.features import mol2alt_sentence, MolSentence, DfVec, sentences2vec, mol2sentence
from mol2vec.helpers import depict_identifier, mol_to_svg, IdentifierTable, plot_2D_vectors

#from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.core.display import HTML
from IPython.display import SVG
IPythonConsole.ipython_useSVG=True
import re




def get_ordered_fingerprint_string_new(SMILE, radius = 2):

    try:
        fingerprints=mol2alt_sentence(Chem.MolFromSmiles(SMILE), radius)
        return ' '.join([str(i) for i in list(fingerprints)])
    
    except Exception as e:
        
        print(e)
        print(SMILE)
        return '' 
    
def get_lenfingerprint(SMILE, radius = 2):

    try:
        return len(set(mol2alt_sentence(Chem.MolFromSmiles(SMILE), radius)))
    
    except Exception as e:
        
        print(e)
        print(SMILE)
        return ''

def get_lenfingerprint2(SMILE, radius = 2):

    try:
        return len(AllChem.GetMorganFingerprint(Chem.MolFromSmiles(SMILE),radius,bitInfo=info).GetNonzeroElements())
    except Exception as e:
        
        print(e)
        print(SMILE)
        return ''
    
    
# explicit representation
def get_binary_representation(SMILE, radius = 2, nBits=2048):
    try:
        array = np.zeros((0, ))
        DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(SMILE),radius, nBits), array)
        return array
    
    except Exception as e:
        
        print(e)
        print(SMILE)
        return ''
    
def get_count_representation(SMILE, radius = 2, nBits=2048):
    try:
        array = np.zeros((0, ))
        DataStructs.ConvertToNumpyArray(AllChem.GetHashedMorganFingerprint(Chem.MolFromSmiles(SMILE),radius, nBits), array)
        return array
    
    except Exception as e:
        
        print(e)
        print(SMILE)
        return ''
    
    
def get_bit_info(SMILE, radius = 2, nBits=2048):   
    try:
        bi = {}
        fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(SMILE), radius, nBits,bitInfo=bi)
        return bi

    except Exception as e:
        print(e)
        print(SMILE)
        return ''    
    

    
    
def get_smile_fragment(SMILE, atomid, radius):   
        mol=Chem.MolFromSmiles(SMILE)
        
        
        if radius >0:
            env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomid)
            atoms=set()
            for bidx in env:
                atoms.add(mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
                atoms.add(mol.GetBondWithIdx(bidx).GetEndAtomIdx())
        
            return Chem.MolFragmentToSmiles(mol,atomsToUse=list(atoms),bondsToUse=env,rootedAtAtom=atomid)
        else:
            return '['+mol.GetAtomWithIdx(atomid).GetSymbol()+']'
        
  

# from rdkit blog


def includeRingMembership(s, n):
    r=';R]'
    d="]"
    return r.join([d.join(s.split(d)[:n]),d.join(s.split(d)[n:])])
 
def includeDegree(s, n, d):
    r=';D'+str(d)+']'
    d="]"
    return r.join([d.join(s.split(d)[:n]),d.join(s.split(d)[n:])])
 
def writePropsToSmiles(mol,smi,order):
    #finalsmi = copy.deepcopy(smi)
    finalsmi = smi
    for i,a in enumerate(order):
        atom = mol.GetAtomWithIdx(a)
        if atom.IsInRing():
            finalsmi = includeRingMembership(finalsmi, i+1)
        finalsmi = includeDegree(finalsmi, i+1, atom.GetDegree())
    return finalsmi

        
def getSubstructSmi(SMILE,atomID,radius):
    mol=Chem.MolFromSmiles(SMILE)
    
    if radius>0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomID)
        atomsToUse=[]
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))
    else:
        atomsToUse = [atomID]
        env=None
    smi = Chem.MolFragmentToSmiles(mol,atomsToUse,bondsToUse=env,allHsExplicit=True, allBondsExplicit=True, rootedAtAtom=atomID)
    order = eval(mol.GetProp("_smilesAtomOutputOrder"))
    smi2 = writePropsToSmiles(mol,smi,order)
    return smi,smi2


def _prepareMol(mol,kekulize):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    return mc
def moltosvg(mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):
    mc = _prepareMol(mol,kekulize)
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc,**kwargs)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return SVG(svg.replace('svg:',''))


# do a depiction where the atom environment is highlighted normally and the central atom
# is highlighted in blue
def getSubstructDepiction(mol,atomID,radius,molSize=(450,200)):
    if radius>0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomID)
        atomsToUse=[]
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))       
    else:
        atomsToUse = [atomID]
        env=None
    return moltosvg(mol,molSize=molSize,highlightAtoms=atomsToUse,highlightAtomColors={atomID:(0.3,0.3,1)})

def depictBit(bitId,SMILE,molSize=(450,200),radius = 2, nBits=2048):
    mol=Chem.MolFromSmiles(SMILE)
    info={}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,radius, nBits,bitInfo=info)
    aid,rad = info[bitId][0]
    return getSubstructDepiction(mol,aid,rad,molSize=molSize)

def depictBit_index(bitId,SMILE,molSize=(450,200),radius = 2, nBits=2048, index=0):
    mol=Chem.MolFromSmiles(SMILE)
    info={}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,radius, nBits,bitInfo=info)
    aid,rad = info[bitId][index]
    return getSubstructDepiction(mol,aid,rad,molSize=molSize)


def NumberAtomsInFragment(SMILE):
    return len(re.sub('[hH]','', "".join(re.findall("[a-zA-Z]+", SMILE))))