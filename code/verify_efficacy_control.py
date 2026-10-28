"""Separately calculate synthetic-control endpoints from saved guide outcomes."""
from pathlib import Path
import hashlib, itertools, json, math
import numpy as np
import pandas as pd

M=Path(__file__).resolve().parents[1]
E=M/'evidence/guide_efficacy_control'
p=json.loads((E/'EFFICACY_CONTROL_PROTOCOL.json').read_text(encoding='utf-8'))
saved=pd.read_csv(E/'synthetic_control_summary.csv')
product_saved=pd.read_csv(E/'all_synthetic_product_selection_contrasts.csv')
arrays=np.load(E/'all_synthetic_guide_outcomes.npz')
checks=[]
for si,sd in enumerate(p['independent_gaussian_noise_sd']):
    for condition in ['shared_efficacy','equal_efficacy','permuted_later_labels']:
        old=arrays[f'{condition}_sd{si}_old'];new=arrays[f'{condition}_sd{si}_new']
        x=old.swapaxes(1,2).reshape(-1,4)
        y=new.swapaxes(1,2).reshape(-1,4)
        xc=x-x.mean(1)[:,None];yc=y-y.mean(1)[:,None]
        denominator=math.sqrt(float(np.einsum('ij,ij->',xc,xc)*np.einsum('ij,ij->',yc,yc)))
        r=float(np.einsum('ij,ij->',xc,yc)/denominator) if denominator>1e-20 else None
        shifts=[]
        plan=[]
        for product in range(len(old)):
            selected=sorted(range(4),key=lambda g:(sum(old[product,g])/3,g))[:2]
            plan.append(selected)
            delta=[]
            for c in range(3):
                all_mean=sum(new[product,:,c])/4
                pairs=[sum(new[product,list(pair),c])/2 for pair in itertools.combinations(range(4),2)]
                assert abs(sum(pairs)/6-all_mean)<1e-12
                delta.append(sum(new[product,selected,c])/2-all_mean)
            shifts.append(sum(delta)/3)
        delta=sum(shifts)/len(shifts)
        ps=product_saved[product_saved.condition.eq(condition)&product_saved.noise_sd.eq(sd)].sort_values('synthetic_product')
        assert ps.synthetic_product.tolist()==list(range(len(old)))
        assert np.array_equal(ps[['selected_guide_1','selected_guide_2']].to_numpy(),plan)
        assert np.allclose(ps.selected_minus_all_mean,shifts,atol=1e-12,rtol=0)
        row=saved[saved.condition.eq(condition)&saved.noise_sd.eq(sd)].iloc[0]
        assert (r is None and pd.isna(row.residual_Pearson_r)) or abs(r-row.residual_Pearson_r)<1e-12
        assert abs(delta-row.selected_minus_all_mean)<1e-12
        checks.append(dict(condition=condition,noise_sd=sd,residual_r=r,mean_shift=delta,
                           all_random_pair_expectations_verified=True))
assert abs(checks[0]['residual_r']-1)<1e-12 and abs(checks[0]['mean_shift']+.3)<1e-12
assert checks[1]['residual_r'] is None and abs(checks[1]['mean_shift'])<1e-12
noiseless=arrays['shared_efficacy_sd0_old']
for g,efficacy in enumerate(p['efficacies']):
    for c,effect in enumerate(p['older_context_effects']):
        assert abs(noiseless[0,g,c]-sum(noiseless[0,:,c])/4
                   -(efficacy-sum(p['efficacies'])/4)*effect)<1e-12
receipt=dict(status='pass',separate_endpoint_implementation=True,conditions_verified=len(checks),
             synthetic_product_selection_contrasts_verified=len(checks)*p['synthetic_products'],
             enumerated_product_context_random_pair_expectations=len(checks)*p['synthetic_products']*3,
             analytical_centring_identity_verified=True,results=checks,
             limitation='Computational verification of an illustrative model; not independent biological validation.')
receipt['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [E/'all_synthetic_guide_outcomes.npz',E/'synthetic_control_summary.csv',E/'all_synthetic_product_selection_contrasts.csv']}
(E/'SEPARATE_CONTROL_VERIFICATION.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in receipt.items() if k!='results'}))
