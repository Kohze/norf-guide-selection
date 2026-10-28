"""Publication masters from frozen evidence; no inference or cohort changes.

180 mm/89 mm masters, Arial, editable vector PDF/SVG and 600 dpi PNG.
All ordinary text >=8 pt; statistics are outside the plotted observations.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch,Rectangle
from matplotlib.text import Text
from matplotlib.markers import MarkerStyle

M=Path(__file__).resolve().parents[1];F=M/'figures';E=M/'evidence/analysis_extensions'
load=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
manifest=load(M/'SUPPLEMENTARY_TABLE_MANIFEST.json')
paths={r['table']:M/r['path'] for r in manifest['tables']}
read=lambda name:pd.read_csv(paths[name],float_precision='round_trip')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
INK='#202C35';BLUE='#0072B2';ORANGE='#D55E00';TEAL='#009E73'
MID='#65717B';LIGHT='#E8EDF0';RULE='#C9D0D5'
CELLS={'A375':(BLUE,'o'),'A549':(ORANGE,'^'),'HEPG2':(TEAL,'s')}
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8.5,
 'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65,
 'xtick.major.size':3,'ytick.major.size':3,'axes.spines.top':False,
 'axes.spines.right':False,'axes.unicode_minus':False,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.default':'regular',
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
 'savefig.facecolor':'white'})
metrics=[]
def drawn_xy(ax):
    arrays=[np.asarray(c.get_offsets(),dtype=float) for c in ax.collections if hasattr(c,'get_offsets') and c.__class__.__name__=='PathCollection']
    return np.concatenate(arrays).tolist() if arrays else []

def figure(height,width=180):return plt.figure(figsize=(width/25.4,height/25.4))
def title(fig,letter,text,y=0.965,x=.02):
    fig.text(x,y,letter,fontweight='bold',fontsize=11,va='top')
    fig.text(x+.045,y,text,fontweight='bold',fontsize=9,va='top')
def txt(ax,x,y,s,**kw):return ax.text(x,y,s,transform=ax.transAxes,**kw)
def zero(ax,x=True,y=False):
    if x:ax.axvline(0,color=MID,lw=.65,ls=(0,(3,3)),zorder=0)
    if y:ax.axhline(0,color=RULE,lw=.6,zorder=0)
def arrow(ax,x,y,xx,yy):ax.add_patch(FancyArrowPatch((x,y),(xx,yy),arrowstyle='-|>',mutation_scale=8,lw=.75,color=MID,shrinkA=0,shrinkB=0))
def save(fig,name,source_tables,**details):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height()
    errors=[];sizes=[]
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():continue
        box=t.get_window_extent(renderer);sizes.append(t.get_fontsize())
        if box.x0<-.5 or box.y0<-.5 or box.x1>w+.5 or box.y1>h+.5:
            errors.append({'text':t.get_text(),'bbox':list(box.extents),'canvas':[w,h]})
    for ext in ['pdf','svg','png']:fig.savefig(F/(name+'.'+ext),dpi=600,facecolor='white')
    metrics.append({'figure':name,'width_mm':float(fig.get_figwidth()*25.4),
      'height_mm':float(fig.get_figheight()*25.4),'minimum_font_pt':min(sizes),
      'text_outside_canvas':errors,'font':'Arial','PNG_dpi':600,
      'source_tables':{n:sha(paths[n]) for n in source_tables},**details})
    plt.close(fig)

# 1: A compact correspondence diagram; branches explicitly use different units.
fig=figure(108);ax=fig.add_axes([.02,.02,.96,.96]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
def section(letter,text,y):
    ax.text(0,y,letter,fontweight='bold',fontsize=11,va='top')
    ax.text(.045,y,text,fontweight='bold',fontsize=9,va='top')
section('a','Complete-product correspondence',.995)
for y,start,color,name in [(.90,.34,BLUE,'Reference product'),(.815,.50,MID,'Alternative start')]:
    ax.text(.025,y,name,va='center',fontsize=8)
    ax.plot([.31,.97],[y,y],color=RULE,lw=.65)
    for left,right in [(start,.65),(.72,.95)]:ax.plot([left,right],[y,y],color=color,lw=5,solid_capstyle='butt')
    ax.plot([start,start],[y-.016,y+.021],color=color,lw=.7)
    ax.text(start,y+.026,'start',ha='center',fontsize=8,color=color)
ax.text(.97,.755,'Illustrative coding chains',ha='right',color=MID,fontsize=8)
section('b','Inventories and constructed relationships',.71)
for x,number,label,note in [(.025,'553 / 2,196','Prensner 2021 / Deutsch 2026','551 older models resolved'),(.40,'274','Constructed relationships','236 from 2021; 268 from 2026'),(.755,'127','Exact-chain/DNA pairs','Each shares ≥1 raw guide')]:
    ax.text(x,.655,number,fontsize=14,fontweight='bold',color=BLUE,va='top')
    ax.text(x,.587,label,fontsize=8,va='top');ax.text(x,.542,note,fontsize=8,color=MID,va='top')
arrow(ax,.31,.607,.37,.607);arrow(ax,.685,.607,.725,.607)
ax.plot([0,1],[.498,.498],color=RULE,lw=.6)
section('c','Distinct cohorts from the 127 exact pairs',.47)
ax.text(.025,.31,'127',fontsize=15,color=BLUE,fontweight='bold',ha='left')
ax.text(.025,.264,'Exact pairs',fontsize=8)
ax.plot([.20,.25,.25],[.315,.315,.411],color=MID,lw=.75)
ax.plot([.25,.25],[.315,.208],color=MID,lw=.75)
for y in [.411,.31,.208]:arrow(ax,.25,y,.30,y)
for y,n,label,note in [(.411,'105','Conditional analysis','471 shared guides'),(.31,'88','Reused / fresh','Separate eligibility gates'),(.208,'16','Additional assay','72 shared guides')]:
    ax.text(.32,y,n,fontsize=11,fontweight='bold',color=INK,va='center')
    ax.text(.40,y,label,fontsize=8.5,va='center')
    ax.text(.40,y-.04,note,fontsize=8,color=MID,va='center')
arrow(ax,.66,.411,.73,.411)
ax.text(.76,.411,'87',fontsize=11,fontweight='bold',color=INK,va='center')
ax.text(.825,.411,'Selection',fontsize=8.5,va='center')
ax.text(.76,.371,'417 shared guides',fontsize=8,color=MID,va='center')
ax.text(.76,.329,'≥4 guides per product',fontsize=8,color=MID,va='center')
ax.plot([0,1],[.145,.145],color=RULE,lw=.6)
section('d','Separate third-publication comparison',.125)
ax.text(.045,.052,'239 sequence/locus matches',fontsize=8.5,va='center')
arrow(ax,.43,.052,.52,.052)
ax.text(.55,.052,'196 products',fontsize=10,fontweight='bold',va='center',color=BLUE)
ax.text(.79,.052,'192 groups',fontsize=8.5,va='center')
save(fig,'Figure1_product_identity_and_cohorts',['S03','S04','S05','S11','S21','S28'])

# 2: All points, shared scales, external statistical annotations and cell shapes.
q=read('S12');a=read('S23')
fig=figure(112);title(fig,'a','Centre each product within each cell and screen')
fig.text(.065,.882,'Guide residual = guide LFC - mean LFC of eligible shared guides',fontsize=8.5)
fig.text(.065,.832,'Equal total weight for each product and fixed cell context',fontsize=8,color=MID)
axes=[fig.add_axes([.105,.19,.365,.49]),fig.add_axes([.61,.19,.365,.49])]
all_values=np.concatenate([q.old_residual,q.new_residual,a.x,a.y]);lo=np.floor((all_values.min()-.12)*2)/2;hi=np.ceil((all_values.max()+.12)*2)/2
for ax,df,x,y,w,letter,heading,counts,stat,old_label in [
 (axes[0],q,'old_residual','new_residual','gene_weight_per_guide_cell','b','Primary screen','105 products; 471 guides','r = 0.643; 95% CI 0.576-0.696','Prensner 2021 primary'),
 (axes[1],a,'x','y','weight','c','Additional assay','16 products; 72 guides','r = 0.911; 95% CI 0.849-0.948','Prensner 2021 additional')]:
    x0=ax.get_position().x0;title(fig,letter,heading,.785,x0-.03)
    fig.text(x0,.722,counts,fontsize=8)
    fig.text(x0,.684,stat,fontsize=8)
    # Deterministically interleave cells so a single context is not always on top.
    plotted=df.iloc[np.random.default_rng(20261026).permutation(len(df))]
    colors=[CELLS[cell][0] for cell in plotted.cell]
    markers=[]
    for cell in plotted.cell:
        shape=MarkerStyle(CELLS[cell][1]);markers.append(shape.get_path().transformed(shape.get_transform()))
    points=ax.scatter(plotted[x],plotted[y],s=7,c=colors,alpha=.65,linewidths=0,zorder=2)
    points.set_paths(markers)
    slope=float((df[w]*df[x]*df[y]).sum()/(df[w]*df[x]*df[x]).sum())
    xx=np.array([float(df[x].min()),float(df[x].max())]);ax.plot(xx,slope*xx,color=INK,lw=.85,zorder=3)
    zero(ax,True,True);ax.set(xlim=(lo,hi),ylim=(lo,hi),xticks=[-2,-1,0,1],yticks=[-2,-1,0,1])
    ax.set_xlabel(old_label+'\nGuide LFC residual (log$_2$ units)',labelpad=6)
axes[0].set_ylabel('Deutsch 2026 guide LFC residual\n(log$_2$ units)',labelpad=6)
handles=[Line2D([],[],marker=marker,color=color,ls='none',markersize=4,label=cell) for cell,(color,marker) in CELLS.items()]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.55,.015),ncol=3,frameon=False,handletextpad=.4,columnspacing=1.8,fontsize=8)
save(fig,'Figure2_conditional_guide_responses',['S12','S23'],main_scatter_points=len(q),additional_scatter_points=len(a),shared_axis_limits=[float(lo),float(hi)],drawn_scatter_xy=[drawn_xy(x) for x in axes])

# 3: Original product-level effect and unchanged post hoc context-exclusion result.
p=read('S27');sens=read('S31');original=load(E/'cycle_93_transferred_guide_selection/GUIDE_SELECTION_RESULTS.json')
fig=figure(107);title(fig,'a','Select in the earlier screen; evaluate in the later screen')
ax=fig.add_axes([.045,.77,.92,.14]);ax.axis('off')
for x,num,text in [(.0,'1','Rank guides by earlier\nequal-cell mean LFC'),(.355,'2','Select two per product\nKeep pair across cells'),(.73,'3','Later selected score minus\nrandom-two expectation')]:
    ax.text(x,.60,num,fontweight='bold',fontsize=11,color=BLUE,va='center')
    ax.text(x+.045,.60,text,fontsize=8,va='center',linespacing=1.45)
arrow(ax,.29,.60,.335,.60);arrow(ax,.66,.60,.71,.60)
title(fig,'b','Original selection',.72,.02);fig.text(.065,.65,'87 products; 79 negative shifts',fontsize=8)
ax=fig.add_axes([.105,.21,.375,.36])
vals=p.newer_selection_shift.to_numpy()
ax.set(xlim=(-1.55,.30),ylim=(-.53,1.18))
# Display-only deterministic swarm in physical point coordinates.
width_pt=fig.get_figwidth()*72*ax.get_position().width
height_pt=fig.get_figheight()*72*ax.get_position().height
xpt=(vals+1.55)/1.85*width_pt;spacing=3.8
placed=[];jitter=np.zeros(len(vals))
for idx in np.argsort(vals,kind='stable'):
    for k in range(200):
        level=0 if k==0 else ((k+1)//2)*(1 if k%2 else -1)
        yp=level*spacing
        if all((xpt[idx]-xp)**2+(yp-y0)**2>=spacing**2-.00001 for xp,y0 in placed):
            placed.append((xpt[idx],yp));jitter[idx]=yp/height_pt*1.71;break
    else:raise RuntimeError('Swarm placement failed')
assert np.max(np.abs(jitter))<.53
ax.scatter(vals,jitter,s=11,color=MID,alpha=.8,linewidths=0)
mu=original['mean_newer_selection_shift'];lo,hi=original['nominal95_cluster_CI']
ax.errorbar(mu,.86,xerr=[[mu-lo],[hi-mu]],fmt='D',color=INK,ms=4,lw=.95,capsize=2.5)
zero(ax);ax.set(xlim=(-1.55,.30),ylim=(-.53,1.18),yticks=[0,.86],yticklabels=['Products','Mean'])
ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
ax.set_xlabel('Later selected-minus-reference LFC\n(log$_2$ units)',labelpad=6)
fig.text(.105,.585,'-0.202 [95% CI -0.246, -0.160]',fontsize=8)
title(fig,'c','Exclude the evaluated cell',.72,.535);fig.text(.58,.65,'Post hoc sensitivity; 95% intervals',fontsize=8)
ax=fig.add_axes([.72,.21,.255,.36])
for y,(cell,(color,marker)) in enumerate(CELLS.items()):
    row=sens[sens.evaluation==cell].iloc[0]
    ax.errorbar(row.mean_shift,y,xerr=[[row.mean_shift-row.nominal95_lower],[row.nominal95_upper-row.mean_shift]],fmt=marker,color=color,ms=4,lw=.95,capsize=2.5)
row=sens[sens.evaluation=='equal_cell_mean'].iloc[0]
ax.errorbar(row.mean_shift,3.35,xerr=[[row.mean_shift-row.nominal95_lower],[row.nominal95_upper-row.mean_shift]],fmt='D',color=INK,ms=4,lw=.95,capsize=2.5)
ax.axhline(2.65,color=RULE,lw=.6);zero(ax)
ax.set(xlim=(-.34,.02),ylim=(3.9,-.65),xticks=[-.3,-.2,-.1,0],yticks=[0,1,2,3.35],yticklabels=['A375','A549','HEPG2','Equal-cell mean'])
ax.tick_params(axis='y',length=0,pad=5);ax.spines['left'].set_visible(False)
ax.set_xlabel('Later selected-minus-reference LFC\n(log$_2$ units)',labelpad=6)
save(fig,'Figure3_transferred_guide_selection',['S27','S31'],product_points=len(p),point_displacement='Deterministic nonoverlapping swarm; x values unchanged',posthoc_intervals=4,drawn_product_shifts=vals.tolist(),swarm_y=jitter.tolist())

# 4: Aligned forest panels, with exact estimates on a separate numerical column.
g=pd.read_csv(M/'figure_source_data/Figure4_product_comparison_estimates.csv',float_precision='round_trip')
fig=figure(94);fig.text(.975,.97,'88 products; nominal 97.5% intervals',fontsize=8,ha='right',va='top',color=MID)
for letter,key,heading,y,limits,ticks,xlabel in [
 ('a','reused_minus_fresh_rank_agreement','Product-level rank agreement',.63,(-.13,.44),[-.1,0,.1,.2,.3,.4],'Reused - fresh Spearman correlation'),
 ('b','fresh_minus_reused_LFC','Mean guide depletion',.20,(-.085,.105),[-.05,0,.05,.1],'Fresh - reused LFC (log$_2$ units)')]:
    title(fig,letter,heading,y+.24)
    ax=fig.add_axes([.335,y,.355,.17]);zero(ax)
    for yy,method,color,marker in [(1,'Observed guide counts',MID,'o'),(0,'Count-balanced (post hoc)',INK,'D')]:
        row=g[(g.endpoint==key)&(g.method==method)].iloc[0]
        ax.errorbar(row.estimate,yy,xerr=[[row.estimate-row.lower],[row.upper-row.estimate]],fmt=marker,color=color,ms=4,lw=.95,capsize=2.5)
        ax.text(-.90,yy,'Observed guide counts' if yy else 'Balanced counts (post hoc)',transform=ax.get_yaxis_transform(),fontsize=8,va='center')
        ax.text(1.10,yy,f'{row.estimate:.3f} [{row.lower:.3f}, {row.upper:.3f}]',transform=ax.get_yaxis_transform(),fontsize=8,va='center')
    ax.set(xlim=limits,ylim=(-.5,1.5),xticks=ticks,yticks=[]);ax.spines['left'].set_visible(False)
    ax.set_xlabel(xlabel,labelpad=5)
    fig.text(.73,y+.20,'Estimate [97.5% CI]',fontsize=8,color=MID)
save(fig,'Figure4_reused_fresh_uncertainty',['S10'],intervals=4,interval_level=.975)

# 5: Product-level evidence with different source baselines explicitly named.
p=read('S29');fig=figure(96)
title(fig,'a','A375 product agreement',.96,.02)
title(fig,'b','Third-publication depletion',.96,.545)
fig.text(.105,.875,'196 products; 192 groups',fontsize=8)
fig.text(.105,.825,'Spearman ρ = 0.462; 95% CI 0.349-0.569',fontsize=8)
fig.text(.625,.875,'0 / 196 at the source strong-depletion cutoff',fontsize=8)
ax=fig.add_axes([.105,.24,.365,.53])
ax.scatter(p.newer_A375_LFC_mean,p.third_A375_LFC_mean2,s=12,color=BLUE,alpha=.65,linewidths=.2,edgecolors='white')
drawn_products=drawn_xy(ax)
zero(ax,True,True);ax.set_xlabel('Deutsch 2026 mean guide LFC\n(plasmid baseline)',labelpad=6)
ax.set_ylabel('Schlesinger 2025 mean product LFC\n(day-4 baseline)',labelpad=6)
ax.set(xticks=[-1,-.5,0,.5],yticks=[-.5,0,.5])
ax=fig.add_axes([.625,.24,.35,.53]);v=p.third_A375_LFC_mean2.to_numpy()
bins=np.linspace(min(-1.10,float(v.min())-.05),max(.65,float(v.max())+.05),18)
counts,_,_=ax.hist(v,bins=bins,color='#739DB5',edgecolor='white',linewidth=.5)
ax.axvline(-1,color=INK,lw=.8,ls=(0,(3,3)));ax.set(xlim=(bins[0],bins[-1]),ylim=(0,max(counts)*1.14),xticks=[-1,-.5,0,.5])
ax.text(-.965,max(counts)*1.07,'LFC = -1',fontsize=8,va='top')
ax.set_xlabel('Schlesinger 2025 mean product LFC\n(day-4 baseline)',labelpad=6);ax.set_ylabel('Products',labelpad=5)
save(fig,'Figure5_third_publication_agreement',['S29'],scatter_points=len(p),histogram_total=int(sum(counts)),histogram_edges=bins.tolist(),histogram_counts=[int(x) for x in counts],drawn_scatter_xy=drawn_products)

# S1: Row-aligned sensitivity estimates; only the primary row has an interval.
r=load(E/'cycle_90_within_product_shared_guide_concordance/WITHIN_PRODUCT_SHARED_GUIDE_CONCORDANCE_RESULTS.json')
names={'min4_shared_guides':'≥4 shared guides','both_screen_splice_safe_min3':'Splice-distance restriction','within_product_cell_guide_ranks':'Within-product guide ranks','descriptive_GC_by_cell_adjusted':'GC adjusted'}
rows=[('Primary estimate',r['primary']['equal_product_cell_residual_Pearson'],105)]
for row in r['secondary']:rows.append((row['cell'] if row['analysis']=='cell_specific' else names[row['analysis']],row['correlation'],row['products']))
fig=figure(92);fig.text(.02,.96,'Sensitivity analysis',fontsize=9,fontweight='bold',va='top')
fig.text(.975,.96,'Only the primary estimate has a nominal 95% interval',fontsize=8,ha='right',va='top',color=MID)
ax=fig.add_axes([.46,.19,.36,.63]);ys=np.arange(len(rows))[::-1]
ax.axhspan(ys[0]-.43,ys[0]+.43,color='#EDF4F8',zorder=-1)
for i,((name,value,n),y) in enumerate(zip(rows,ys)):
    if i==0:
        lo,hi=r['primary']['cluster_bootstrap_CI95'];ax.errorbar(value,y,xerr=[[value-lo],[hi-value]],fmt='D',color=INK,ms=4,lw=.95,capsize=2.5)
    else:
        color,marker=CELLS.get(name,(MID,'o'));ax.plot(value,y,marker,color=color,ms=4)
    ax.text(-1.22,y,name,transform=ax.get_yaxis_transform(),fontsize=8,va='center',fontweight='bold' if i==0 else 'normal')
    ax.text(-.24,y,str(n),transform=ax.get_yaxis_transform(),fontsize=8,va='center',ha='right')
    ax.text(1.28,y,f'{value:.3f}',transform=ax.get_yaxis_transform(),fontsize=8,va='center',ha='right')
ax.set(xlim=(.40,.85),ylim=(-.6,7.6),yticks=[],xticks=[.4,.5,.6,.7,.8]);ax.spines['left'].set_visible(False)
ax.grid(axis='x',color=LIGHT,lw=.6);ax.set_xlabel('Weighted residual correlation',labelpad=6)
fig.text(.02,.85,'Analysis',fontsize=8,color=MID);fig.text(.374,.85,'Products',fontsize=8,color=MID,ha='right');fig.text(.92,.85,'r',fontsize=8,color=MID,ha='right')
save(fig,'FigureS1_conditional_sensitivities',['S13'],sensitivity_points=len(rows),inferential_intervals=1)

# S2: Single-column descriptive comparison, with a genuine zero baseline.
t=load(E/'cycle_94_independent_product_replication/PRODUCT_CONCORDANCE_RESULTS.json')
fig=figure(77,width=89);ax=fig.add_axes([.15,.22,.78,.55])
values=[t['top20percent']['overlap'],t['top20percent']['uniform_random_expectation']]
ax.bar([0,1],values,color=[BLUE,'#AEBCC6'],width=.52,linewidth=0)
for x,value in enumerate(values):ax.text(x,value+.6,f'{value:.2f}'.rstrip('0').rstrip('.'),ha='center',fontsize=9,fontweight='bold')
ax.set(ylim=(0,22),yticks=[0,5,10,15,20],xticks=[0,1],xticklabels=['Observed','Uniform independent\nexpectation'])
ax.set_ylabel('Products shared by top-40 sets',labelpad=5)
fig.text(.15,.94,'196 products; 40 selected per study',fontsize=8,va='top')
fig.text(.15,.87,'Descriptive comparison',fontsize=8,color=MID,va='top')
save(fig,'FigureS2_third_publication_priority_overlap',['S29'],width_role='single_column',bar_values=values,axis_starts_at_zero=True)

receipt={'status':'rendered_pending_visual_review','new_outcome_analysis':False,
 'new_intervals_generated':False,'masters':metrics,'ordinary_font_min_pt':8,
 'formats':['editable SVG','vector PDF with TrueType fonts','600 dpi PNG'],
 'all_scatter_observations_retained':True}
(M/'FIGURE_BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
bad=[{'figure':m['figure'],'errors':m['text_outside_canvas']} for m in metrics if m['text_outside_canvas']]
print(json.dumps({'figures':len(metrics),'outside_canvas':bad,'min_font_sizes':{m['figure']:m['minimum_font_pt'] for m in metrics}}))
