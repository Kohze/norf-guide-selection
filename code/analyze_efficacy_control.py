"""Fixed synthetic positive control, not an empirical efficacy estimate.

Protocol: ../evidence/guide_efficacy_control/EFFICACY_CONTROL_PROTOCOL.json.
No empirical guide outcomes are loaded by this program.
"""
from pathlib import Path
import hashlib, itertools, json
import numpy as np
import pandas as pd

M = Path(__file__).resolve().parents[1]
E = M / 'evidence/guide_efficacy_control'
E.mkdir(exist_ok=True)
protocol = json.loads((E/'EFFICACY_CONTROL_PROTOCOL.json').read_text(encoding='utf-8'))
N = protocol['synthetic_products']
C = protocol['fixed_contexts']
G = protocol['guides_per_synthetic_product']
a = np.asarray(protocol['efficacies'])
theta_old = np.asarray(protocol['older_context_effects'])
theta_new = np.asarray(protocol['later_context_effects'])
rng = np.random.default_rng(protocol['seed'])
conditions = ['shared_efficacy', 'equal_efficacy', 'permuted_later_labels']
aggregates, products, arrays = [], [], {}
for si, sd in enumerate(protocol['independent_gaussian_noise_sd']):
    noise_old = rng.normal(0, sd, (N,G,C))
    noise_new = rng.normal(0, sd, (N,G,C))
    heterogeneous_old = a[None,:,None]*theta_old[None,None,:] + noise_old
    heterogeneous_new = a[None,:,None]*theta_new[None,None,:] + noise_new
    permutations = np.asarray([rng.permutation(G) for _ in range(N)])
    for condition in conditions:
        if condition == 'equal_efficacy':
            old = np.mean(a)*theta_old[None,None,:] + noise_old
            new = np.mean(a)*theta_new[None,None,:] + noise_new
        else:
            old = heterogeneous_old.copy()
            new = heterogeneous_new.copy()
            if condition == 'permuted_later_labels':
                new = np.take_along_axis(new, permutations[:,:,None], axis=1)
        ex = old - old.mean(axis=1,keepdims=True)
        ey = new - new.mean(axis=1,keepdims=True)
        xx, yy, xy = np.mean(ex*ex), np.mean(ey*ey), np.mean(ex*ey)
        r = float(xy/np.sqrt(xx*yy)) if xx > 1e-25 and yy > 1e-25 else None
        # Stable sorting resolves ties by the synthetic guide index.
        selected = np.argsort(old.mean(axis=2),axis=1,kind='stable')[:,:2]
        selected_new = np.take_along_axis(new, selected[:,:,None], axis=1).mean(axis=1)
        per_product = (selected_new-new.mean(axis=1)).mean(axis=1)
        shift = float(per_product.mean())
        aggregates.append(dict(condition=condition,noise_sd=sd,synthetic_products=N,
                               guides_per_product=G,fixed_contexts=C,
                               residual_Pearson_r=r,selected_minus_all_mean=shift))
        for p in range(N):
            products.append(dict(condition=condition,noise_sd=sd,synthetic_product=p,
                                 selected_guide_1=int(selected[p,0]),selected_guide_2=int(selected[p,1]),
                                 selected_minus_all_mean=float(per_product[p])))
        key=f'{condition}_sd{si}'
        arrays[key+'_old']=old
        arrays[key+'_new']=new
        # Exact expectation over all equally likely two-guide plans.
        all_pairs = np.stack([new[:,list(pair)].mean(axis=1)
                              for pair in itertools.combinations(range(G),2)])
        assert np.allclose(all_pairs.mean(axis=0),new.mean(axis=1),atol=1e-14)
        if sd == 0 and condition == 'shared_efficacy':
            assert np.allclose(ex,(a-a.mean())[None,:,None]*theta_old[None,None,:])
            assert abs(r-1)<1e-12
            assert np.allclose(per_product,-0.3)
        if sd == 0 and condition == 'equal_efficacy':
            assert r is None and np.allclose(per_product,0)
np.savez_compressed(E/'all_synthetic_guide_outcomes.npz',**arrays)
summary = pd.DataFrame(aggregates)
summary.to_csv(E/'synthetic_control_summary.csv',index=False)
pd.DataFrame(products).to_csv(E/'all_synthetic_product_selection_contrasts.csv',index=False)
receipt = dict(kind='fixed_parameter_synthetic_positive_control',parameters_fitted_to_empirical_data=False,
               empirical_guide_outcomes_read=False,empirical_mechanism_identified=False,
               empirical_variance_or_covariance_explained_estimated=False,
               synthetic_products_per_condition=N,conditions=len(conditions),noise_levels=4,
               analytical_noiseless_shared_r=1.0,analytical_noiseless_shared_shift=-0.3,
               interpretation='Shared on-target efficacy acting on true dependencies can generate both endpoints. This demonstrates compatibility; it does not estimate whether efficacy explains the observed data.',
               results=aggregates,
               output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in E.iterdir() if p.suffix in ['.csv','.npz']})
(E/'EFFICACY_CONTROL_RESULTS.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(summary.to_string(index=False))

# Publication artwork: editable text, vector PDF and 600-dpi PNG.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
                     'xtick.labelsize':8,'ytick.labelsize':8,
                     'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'stix',
                     'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(180/25.4,83/25.4))
fig.subplots_adjust(left=.085,right=.985,bottom=.29,top=.89,wspace=.34)
spec=[('shared_efficacy','Shared efficacy','#0072B2','o','-'),
      ('equal_efficacy','Equal efficacy','#666666','s','--'),
      ('permuted_later_labels','Permuted later labels','#D55E00','^',':')]
for key,label,color,marker,ls in spec:
    d=summary[summary.condition.eq(key)]
    for ax,col in zip(axes,['residual_Pearson_r','selected_minus_all_mean']):
        ax.plot(d.noise_sd,d[col],label=label,color=color,marker=marker,ms=4,
                linewidth=1.1,linestyle=ls,clip_on=False)
for ax,title in zip(axes,['Centred guide responses','Transferred selection']):
    ax.set_title(title,fontsize=9,loc='left',pad=11)
    ax.set_xlabel('Independent noise SD (log$_2$ units)',labelpad=6)
    ax.set_xticks([0,.25,.5,1]);ax.set_xlim(-.035,1.035)
    ax.axhline(0,color='#bfc5cb',linewidth=.7,zorder=0)
    ax.tick_params(length=3,width=.6);ax.spines['left'].set_linewidth(.6);ax.spines['bottom'].set_linewidth(.6)
axes[0].set_ylim(-.08,1.06);axes[0].set_yticks([0,.25,.5,.75,1])
axes[0].set_ylabel('Residual Pearson correlation')
axes[1].set_ylim(-.325,.025);axes[1].set_yticks([-.3,-.2,-.1,0])
axes[1].set_ylabel('Selected minus all-guide mean (log$_2$)')
axes[0].text(-.15,1.07,'a',transform=axes[0].transAxes,fontweight='bold',fontsize=10)
axes[1].text(-.15,1.07,'b',transform=axes[1].transAxes,fontweight='bold',fontsize=10)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.065),
           ncol=3,frameon=False,columnspacing=1.7,handlelength=2)
fig.text(.5,.01,'Fixed synthetic positive control; no empirical parameters fitted',ha='center',fontsize=8,color='#4c5660')
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
for t in fig.findobj(matplotlib.text.Text):
    if not t.get_visible() or not t.get_text(): continue
    bb=t.get_window_extent(renderer)
    if t in fig.texts or t in fig.legends or t.get_text() in ['a','b']:
        assert bb.x0>=0 and bb.y0>=0 and bb.x1<=fig.bbox.width and bb.y1<=fig.bbox.height,t.get_text()
out=M/'figures/FigureS3_efficacy_positive_control'
out.parent.mkdir(exist_ok=True)
for suffix in ['.pdf','.svg','.png']:
    fig.savefig(out.with_suffix(suffix),dpi=600)
plt.close(fig)
