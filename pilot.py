"""AAV9 local, reference-site-informed pilot. Not blind docking.
Run from project root: python src/pilot.py prepare|dock|analyze
"""
import argparse, json, pathlib, time, platform, importlib.metadata, csv
import numpy as np
import gemmi
from scipy.spatial.distance import cdist
R=pathlib.Path(__file__).resolve().parents[1]
D=R/'data/processed'; O=R/'results'; D.mkdir(exist_ok=True); O.mkdir(exist_ok=True)
def xyz(res): return np.array([list(a.pos) for a in res if not a.element.is_hydrogen])
def protein(s): return {c.name:c for c in s[0] if any(r.het_flag=='A' for r in c)}
def atomdict(chains):
 return {(cn,r.seqid.num,r.name,a.name):np.array(list(a.pos)) for cn,c in chains.items() for r in c if r.het_flag=='A' for a in r if not a.element.is_hydrogen}
def fit(x,y):
 xc=x.mean(0);yc=y.mean(0);u,_,vt=np.linalg.svd((x-xc).T@(y-yc));rot=u@vt
 if np.linalg.det(rot)<0:u[:,-1]*=-1;rot=u@vt
 return rot,yc-xc@rot

def prepare():
 h=gemmi.read_structure(str(R/'data/raw/7MTZ.cif'));a=gemmi.read_structure(str(R/'data/raw/7MT0.cif'))
 lig=next(r for c in h[0] for r in c if r.name=='GAL'); ref=xyz(lig);center=ref.mean(0)
 hp=protein(h);ap=protein(a);hd=atomdict(hp);ad=atomdict(ap)
 # Entire subunits touching a 20-A ligand shell; retains neighbors around a 20-A box.
 selected=[cn for cn,c in hp.items() if cdist(np.vstack([xyz(r) for r in c if r.het_flag=='A']),ref).min()<20]
 keys=[k for k in hd if k[0] in selected and k[3]=='CA' and k in ad]
 rot,trans=fit(np.array([ad[k] for k in keys]),np.array([hd[k] for k in keys]))
 rmsd=float(np.sqrt(np.mean(np.sum((np.array([ad[k] for k in keys])@rot+trans-np.array([hd[k] for k in keys]))**2,axis=1))))
 # Preserve author numbering, remap selected chains to short IDs for legacy PDB tools.
 chainmap={cn:chr(65+i) for i,cn in enumerate(selected)}
 for name,pc in [('holo',hp),('apo',ap)]:
  out=gemmi.Structure();out.add_model(gemmi.Model('1'))
  for cn in selected:
   c=pc[cn].clone();c.name=chainmap[cn]
   for i in reversed(range(len(c))):
    if c[i].het_flag!='A':del c[i]
   if name=='apo':
    for res in c:
     for atom in res:
      p=np.array(list(atom.pos))@rot+trans;atom.pos=gemmi.Position(*p)
   out[0].add_chain(c)
  out.write_pdb(str(D/f'{name}_receptor.pdb'))
 contacts=[]
 for k,p in hd.items():
  ds=np.linalg.norm(ref-p,axis=1)
  for i in np.where(ds<4.0)[0]: contacts.append({'chain':k[0],'resnum':k[1],'resname':k[2],'atom':k[3],'ligand_atom':lig[i].name,'distance_A':float(ds[i])})
 # Build ligand from CCD explicit connectivity; never infer bonds from pose proximity.
 from rdkit import Chem
 from rdkit.Chem import AllChem
 from meeko import MoleculePreparation,PDBQTWriterLegacy
 b=gemmi.cif.read_file(str(R/'data/raw/GAL.cif')).sole_block();m=Chem.RWMol();ids={}
 for row in b.find('_chem_comp_atom.',['atom_id','type_symbol','charge']):
  atom=Chem.Atom(row[1]);atom.SetFormalCharge(int(row[2]));atom.SetProp('ccd_name',row[0]);ids[row[0]]=m.AddAtom(atom)
 for row in b.find('_chem_comp_bond.',['atom_id_1','atom_id_2','value_order']):m.AddBond(ids[row[0]],ids[row[1]],{'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE}[row[2]])
 mol=m.GetMol();Chem.SanitizeMol(mol)
 # Set ideal CCD coordinates (including H), assign stereochemistry there, then embed independently.
 conf=Chem.Conformer(mol.GetNumAtoms())
 for row in b.find('_chem_comp_atom.',['atom_id','pdbx_model_Cartn_x_ideal','pdbx_model_Cartn_y_ideal','pdbx_model_Cartn_z_ideal']):conf.SetAtomPosition(ids[row[0]],tuple(map(float,list(row)[1:])))
 mol.AddConformer(conf);Chem.AssignAtomChiralTagsFromStructure(mol);Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
 initial_smiles=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)
 mol.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=20261024
 if AllChem.EmbedMolecule(mol,params)!=0:raise RuntimeError('ligand embedding failed')
 AllChem.MMFFOptimizeMolecule(mol)
 conf=mol.GetConformer();pos=conf.GetPositions();pos+=center-pos.mean(0)
 for i,p in enumerate(pos):conf.SetAtomPosition(i,p)
 # Names survive Meeko PDBQT output via PDB residue info.
 for at in mol.GetAtoms():
  info=Chem.AtomPDBResidueInfo();info.SetName(at.GetProp('ccd_name').rjust(4));info.SetResidueName('GAL');info.SetResidueNumber(1);at.SetMonomerInfo(info)
 setup=MoleculePreparation().prepare(mol)[0];qt,ok,err=PDBQTWriterLegacy.write_string(setup)
 if not ok:raise RuntimeError(err)
 (D/'galactose.pdbqt').write_text(qt)
 w=Chem.SDWriter(str(D/'galactose.sdf'));w.write(mol);w.close()
 (D/'reference.json').write_text(json.dumps({'names':[at.name for at in lig],'coordinates':ref.tolist()},indent=2))
 result={'status':'structural_preparation_complete','selected_author_chains':selected,'pdb_chain_mapping':chainmap,'box_center':center.tolist(),'box_size':[20,20,20],'selection_shell_A':20,'apo_to_holo_CA_rmsd_A':rmsd,'alignment_CA_count':len(keys),'rotation':rot.tolist(),'translation':trans.tolist(),'galactose_isomeric_smiles':initial_smiles,'contacts_4A':contacts,'contact_cutoff_A':4.0,'reference_ligand':'7MTZ author chain A GAL 801','information_regime':'true-site-local; not blind','versions':{x:importlib.metadata.version(x) for x in ['gemmi','rdkit','vina','meeko','numpy','scipy']},'platform':platform.platform()}
 (O/'preparation.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='contacts_4A'},indent=2))

def dock():
 from vina import Vina
 meta=json.loads((O/'preparation.json').read_text());runs=[]
 for receptor in ['holo','apo']:
  if not (D/f'{receptor}_receptor.pdbqt').exists():
   (O/f'{receptor}_blocked.json').write_text(json.dumps({'status':'blocked','reason':'prepared receptor PDBQT missing; inspect preparation log'}));continue
  for seed in [104729,130363,155921]:
   start=time.perf_counter();v=Vina(sf_name='vina',cpu=2,seed=seed,verbosity=0)
   v.set_receptor(str(D/f'{receptor}_receptor.pdbqt'));v.set_ligand_from_file(str(D/'galactose.pdbqt'))
   v.compute_vina_maps(center=meta['box_center'],box_size=meta['box_size']);v.dock(exhaustiveness=16,n_poses=20,min_rmsd=0.5)
   path=O/f'{receptor}_{seed}.pdbqt';v.write_poses(str(path),n_poses=20,energy_range=10,overwrite=True)
   runs.append({'receptor':receptor,'seed':seed,'wall_seconds':time.perf_counter()-start,'exhaustiveness':16,'cpu':2,'requested_poses':20,'pose_file':path.name});(O/'runs.json').write_text(json.dumps(runs,indent=2));print(runs[-1],flush=True)

def analyze():
 ref=json.loads((D/'reference.json').read_text());ref=dict(zip(ref['names'],ref['coordinates']));names=list(ref);target=np.array([ref[n] for n in names]);rows=[]
 for run in json.loads((O/'runs.json').read_text()):
  receptor=gemmi.read_structure(str(D/f"{run['receptor']}_receptor.pdb"));pd=atomdict(protein(receptor));pk=list(pd);px=np.array(list(pd.values()))
  holo=atomdict(protein(gemmi.read_structure(str(D/'holo_receptor.pdb')))); hx=np.array([holo[k] for k in pk]);native=set(map(tuple,np.argwhere(cdist(target,hx)<4)))
  for rank,block in enumerate((O/run['pose_file']).read_text().split('MODEL')[1:],1):
   atomlines=[l for l in block.splitlines() if l.startswith(('ATOM','HETATM'))];coords={l[12:16].strip():[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in atomlines if l[12:16].strip() in ref}
   if set(coords)!=set(ref):raise ValueError('Atom identity mismatch; refusing RMSD')
   x=np.array([coords[n] for n in names]);rmsd=float(np.sqrt(np.mean(np.sum((x-target)**2,axis=1))));pred=set(map(tuple,np.argwhere(cdist(x,px)<4)))
   # Each receptor compared with the same reference ligand in its aligned frame.
   score=float(next(l for l in block.splitlines() if l.startswith('REMARK VINA RESULT:')).split()[3])
   rows.append({**run,'rank':rank,'vina_score':score,'rmsd_A':rmsd,'contact_recall':len(native&pred)/len(native) if native else None,'contact_precision':len(native&pred)/len(pred) if pred else 0,'success_2A':rmsd<=2,'success_3A':rmsd<=3})
 with (O/'pose_metrics.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 summary=[]
 for run in json.loads((O/'runs.json').read_text()):
  rr=[x for x in rows if x['pose_file']==run['pose_file']];best=min(rr,key=lambda x:x['rmsd_A'])
  summary.append({**run,'returned_poses':len(rr),'top1_rmsd_A':rr[0]['rmsd_A'],'top5_min_rmsd_A':min(x['rmsd_A'] for x in rr[:5]),'best_returned_rmsd_A':best['rmsd_A'],'best_returned_rank':best['rank'],'top1_score':rr[0]['vina_score']})
 (O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','dock','analyze']);globals()[p.parse_args().stage]()
