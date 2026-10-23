"""Coordinate geometry; all distances are Angstrom, all torsions are degrees.

Alpha-like torsions and O--N distances are descriptive proxies, never DSSP or
evidence of hydrogen-bond occupancy. Direct comparison is used in frozen runs.
"""
from __future__ import annotations
import math
import re
import numpy as np

BACKBONE = frozenset(('N', 'CA', 'C', 'O', 'OXT'))


def coordinates(text: str) -> dict:
    result = {}
    for line in text.splitlines():
        if not line.startswith('ATOM  '):
            continue
        name = line[12:16].strip()
        element = (line[76:78].strip() if len(line) >= 78 else '') or re.sub(r'^[0-9]+', '', name)[0]
        if element.upper() in ('H', 'D'):
            continue
        key = (line[21], int(line[22:26]), name)
        if key in result:
            raise ValueError('Duplicate atom identity')
        xyz = tuple(float(line[k:k+8]) for k in (30,38,46))
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError('Nonfinite coordinate')
        result[key] = xyz
    return result


def select(coords: dict, chain: str, backbone=False, start=None, end=None) -> dict:
    return {k:v for k,v in coords.items() if k[0] == chain and (not backbone or k[2] in BACKBONE) and (start is None or k[1] >= start) and (end is None or k[1] <= end)}


def displacement(reference: dict, mobile: dict) -> dict:
    if not reference or reference.keys() != mobile.keys():
        raise ValueError('Missing/mismatched atom identities')
    delta=np.asarray([mobile[k] for k in reference])-np.asarray(list(reference.values()))
    norm=np.linalg.norm(delta,axis=1)
    return dict(rmsd_A=float(np.sqrt(np.mean(norm**2))),max_A=float(norm.max()),n_atoms=len(norm))


def align_on_receptor(reference: dict, mobile: dict) -> dict:
    keys=sorted(k for k in reference if k[0]=='A' and k[2]=='CA')
    if len(keys)<3 or any(k not in mobile for k in keys):
        raise ValueError('Insufficient receptor CA correspondence')
    target=np.asarray([reference[k] for k in keys]); moving=np.asarray([mobile[k] for k in keys])
    tc,mc=target.mean(0),moving.mean(0)
    u,_,vt=np.linalg.svd((moving-mc).T@(target-tc))
    correction=np.eye(3); correction[-1,-1]=np.linalg.det(u@vt)
    rotation=u@correction@vt
    return {k:tuple((np.asarray(v)-mc)@rotation+tc) for k,v in mobile.items()}


def dihedral(p0,p1,p2,p3) -> float:
    p0,p1,p2,p3=map(lambda p:np.asarray(p,dtype=float),(p0,p1,p2,p3))
    b0=-(p1-p0); b1=p2-p1; b2=p3-p2
    norm=np.linalg.norm(b1)
    if norm==0: raise ValueError('Degenerate torsion')
    b1=b1/norm; v=b0-np.dot(b0,b1)*b1; w=b2-np.dot(b2,b1)*b1
    if np.linalg.norm(v)<1e-12 or np.linalg.norm(w)<1e-12: raise ValueError('Degenerate torsion')
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w))))


def torsions(coords: dict, sequence: str) -> list[dict]:
    rows=[]
    for i,aa in enumerate(sequence,1):
        phi=dihedral(coords[('B',i-1,'C')],coords[('B',i,'N')],coords[('B',i,'CA')],coords[('B',i,'C')]) if i>1 else None
        psi=dihedral(coords[('B',i,'N')],coords[('B',i,'CA')],coords[('B',i,'C')],coords[('B',i+1,'N')]) if i<len(sequence) else None
        rows.append(dict(residue=i,amino_acid=aa,phi_deg=phi,psi_deg=psi,phi_status='defined' if phi is not None else 'terminal_missing',psi_status='defined' if psi is not None else 'terminal_missing',interior_2_11=2<=i<=11,local_8_11=8<=i<=11,alpha_like=phi is not None and psi is not None and -100<=phi<=-30 and -80<=psi<=-5))
    return rows


def hbond_proxies(coords: dict, sequence: str) -> list[dict]:
    rows=[]
    for i in range(1,len(sequence)-3):
        donor=i+4
        distance=float(np.linalg.norm(np.asarray(coords[('B',i,'O')])-np.asarray(coords[('B',donor,'N')])))
        eligible=sequence[donor-1]!='P'
        rows.append(dict(acceptor_residue=i,donor_residue=donor,distance_A=distance,eligible_donor=eligible,proxy_present=eligible and distance<=3.5,description='heavy_atom_distance_proxy_not_Hbond_occupancy'))
    return rows


def geometry(coords: dict, sequence: str) -> dict:
    t=torsions(coords,sequence); h=hbond_proxies(coords,sequence)
    return dict(alpha_like_2_11=sum(x['alpha_like'] for x in t if x['interior_2_11']),alpha_like_8_11=sum(x['alpha_like'] for x in t if x['local_8_11']),hbond_distance_proxy_count=sum(x['proxy_present'] for x in h))


def movement(reference: dict, mobile: dict) -> dict:
    result={}
    for prefix,chain,bb,start,end in [('A_backbone','A',True,None,None),('B_backbone','B',True,None,None),('B_all','B',False,None,None),('B8_12_backbone','B',True,8,12),('B8_12_all','B',False,8,12)]:
        metric=displacement(select(reference,chain,bb,start,end),select(mobile,chain,bb,start,end))
        result.update({prefix+'_'+k:v for k,v in metric.items()})
    return result
