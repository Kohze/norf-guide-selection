from pathlib import Path
import csv,json,gzip,math,statistics,hashlib,datetime
import numpy as np
O=Path(__file__).resolve().parent;N=O.parent
C66=N/'cycle_66_independent_nORF_fitness_benchmark';C83=N/'cycle_83_human_nORF_guide_splice_liability';C35=N/'cycle_35_guide_replicates_and_medical_provenance';C34=N/'cycle_34_endogenous_peptide_and_ORF_provenance'
def rows(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def truth(x):return str(x).lower()=='true'
checks=0;maxerr=0
def chk(x): 
 global checks
 assert x;checks+=1
def eq(a,b,tol=2e-11):
 global checks,maxerr
 a=float(a);b=float(b)
 if math.isnan(a) or math.isnan(b):chk(math.isnan(a) and math.isnan(b));return
 err=abs(a-b);maxerr=max(maxerr,err);chk(err<=tol*(1+abs(b)))
def mean(v):return math.fsum(v)/len(v)
def ranks(v):
 order=sorted(range(len(v)),key=lambda i:v[i]);out=[0.]*len(v);j=0
 while j<len(v):
  k=j+1
  while k<len(v) and v[order[k]]==v[order[j]]:k+=1
  for q in range(j,k):out[order[q]]=(j+1+k)/2
  j=k
 return np.array(out)
def rho(a,b):
 if len(a)<3 or len(set(a))<2 or len(set(b))<2:return float('nan')
 x=ranks(a);y=ranks(b);x=x-x.mean();y=y-y.mean();return float(np.dot(x,y)/math.sqrt(np.dot(x,x)*np.dot(y,y)))
def stats(A,R,F,retain=(0,1,2)):
 return (mean((F[:,retain]-R[:,retain]).flatten()),mean([rho(A[:,c],R[:,c])-rho(A[:,c],F[:,c]) for c in retain]))
pins=json.loads((O/'PRE_OUTCOME_JOIN_INPUT_AND_SOURCE_PINS.json').read_text())['pins']
for p in pins:
 b=Path(p['path']).read_bytes();chk(len(b)==p['bytes']);chk(hashlib.sha256(b).hexdigest()==p['sha256'])
cells=['A375','A549','HEPG2']
og=rows(O/'all1010_older_exact_product_guide_eligibility_before_outcome_join.csv');ng=rows(O/'all720_newer_exact_product_guide_roles_and_eligibility_before_outcome_join.csv')
oldgeom={r['sgrna']:r for r in rows(C83/'all_original2021_library_guide_splice_and_CDS_geometry_before_outcome_join.csv')}
newgeom={r['sgrna']:r for r in rows(C83/'all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv')}
oldlib=rows(C83/'original2021_complete_guide_library_for_genome_scan.csv');oldall={r['sgrna'] for r in oldlib};oldassigned={}
for r in oldlib:oldassigned.setdefault(r['gene'],set()).add(r['sgrna'])
rawcounts={r['sgrna']:r for r in rows(C35/'all_original_table13_B_guide_measurements.csv')}
oldout={(r['Construct Barcode'],r['Cell Line']):r for r in rows(C66/'all_original_primary_screen_guide_LFC_rows.csv')}
newout={(r['sgrna'],r['cell']):r for r in rows(C83/'all215704_original2026_guide_cell_primary_log_abundance_changes.csv')}
for group,geo,outcomes in [(og,oldgeom,oldout),(ng,newgeom,newout)]:
 for q in group:
  src=geo[q['sgrna']]
  for k in src:chk(q[k]==src[k])
  unique=truth(src['unique_primary_NGG_site']) and truth(src['cut_inside_assigned_source_ORF']) and int(src['source_distinct_target_assignments_for_guide'])==1
  present=all((q['sgrna'],c) in outcomes for c in cells);chk(truth(q['source_present_all3sharedcells'])==present)
  coverage=True if group is og else all(float(rawcounts[q['sgrna']]['pDNA_'+str(k)])>=10 for k in (1,2,3))
  chk(truth(q['unique_only_eligible'])==(unique and present and coverage))
  chk(truth(q['primary_eligible'])==(truth(src['primary_geometry_eligible']) and present and coverage))
  if group is ng:
   role='reused' if q['sgrna'] in oldassigned[q['old_source_name']] else ('other_old_library_reuse' if q['sgrna'] in oldall else 'fresh')
   chk(q['guide_group']==role);chk(truth(q['plasmid_coverage_pass'])==coverage)
oldv={(q,c):float(r['LFC_D21_relative_to_pDNA']) for (q,c),r in oldout.items() if c in cells}
# Independent scalar logarithm source calculation; archived values retained only after this comparison.
newv={};scalarv={}
for q in {r['sgrna'] for r in ng}:
 row=rawcounts[q]
 for c in cells:
  late=[float(row[c+'-311Cas9_Rep '+k]) for k in 'ABC'];pl=[float(row['pDNA_'+str(k)]) for k in (1,2,3)]
  for pc in (.5,1.,5.):
   v=(math.fsum(math.log2(x+pc) for x in late)-math.fsum(math.log2(x+pc) for x in pl))/3
   scalarv[(q,c,pc)]=v
   nv=float(np.log2(np.array(late)+pc).mean()-np.log2(np.array(pl)+pc).mean())
   eq(v,nv)
   newv[(q,c,pc)]=nv
  eq(scalarv[(q,c,1.)],newout[(q,c)]['LFC'])
for r in rows(O/'all3030_older_exact_product_source_guide_cell_outcomes.csv'):eq(r['LFC'],oldv[(r['sgrna'],r['cell'])])
for r in rows(O/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv'):eq(r['LFC'],newout[(r['sgrna'],r['cell'])]['LFC'])
for r in rows(O/'all4320_newer_guide_cell_alternative_pseudocount_effects.csv'):eq(r['LFC'],scalarv[(r['sgrna'],r['cell'],float(r['pseudocount']))])
den=rows(O/'all127_exact_product_cohort_eligibility_and_clusters_before_outcome_join.csv');dm={r['pair_ID']:r for r in den}
oldby={i:[r for r in og if r['pair_ID']==i] for i in dm};newby={i:[r for r in ng if r['pair_ID']==i] for i in dm}
# Independent connected-component check of both supplied dependency partitions.
oc={r['source_Excel_row']:r['source_locus_cluster'] for r in rows(C66/'all553_analysis_source_labels_context_and_locus_groups.csv')}
nc={r['orf_targeted']:r['dependence_union_cluster'] for r in rows(N/'cycle_70_ORF_resolved_context_and_host_dependency'/'all17568_source_resolved_context_host_fitness_and_eligibility.csv')}
tokens={i:{('old',oc[r['older_source_Excel_row']]),('new',nc[r['newer_source_ORF_ID']])} for i,r in dm.items()}
for i in dm:
 reached={i};changed=True
 while changed:
  changed=False
  for j in dm:
   if j not in reached and any(tokens[j]&tokens[k] for k in reached):reached.add(j);changed=True
 for j in dm:chk((j in reached)==(dm[i]['dependence_cluster']==dm[j]['dependence_cluster']))
expr={(r['orf_targeted'],r['cell_line']):r for r in rows(C34/'table13_new_eight_cell_fitness.csv')}
def numeric(x):
 try:return float(x)
 except:return float('nan')
for i,d in dm.items():
 for tag in ('primary','unique_only'):
  flag=tag+'_eligible';a=[q for q in oldby[i] if truth(q[flag])];b=[q for q in newby[i] if truth(q[flag])]
  nr=sum(q['guide_group']=='reused' for q in b);nf=sum(q['guide_group']=='fresh' for q in b)
  for key,v in [('older_guides',len(a)),('reused_guides',nr),('fresh_guides',nf)]:chk(int(d[tag+'_'+key])==v)
  chk(truth(d[flag])==(len(a)>0 and nr>0 and nf>0))
 er=[expr[(d['newer_source_ORF_ID'],c)] for c in cells]
 known=all(math.isfinite(numeric(e['rnaseq_tpm'])) and math.isfinite(numeric(e['riboseq_tpm'])) for e in er)
 qualified=known and all(numeric(e['rnaseq_tpm'])>=10 and numeric(e['riboseq_tpm'])>=5 for e in er)
 chk(truth(d['expression_known_all3cells'])==known);chk(truth(d['expression_qualified_all3cells'])==qualified)
 chk(truth(d['expression_qualified_cohort_eligible'])==(qualified and truth(d['primary_eligible'])))
 chk(truth(d['source_AA_clean_eligible'])==(truth(d['primary_eligible']) and not truth(d['prior_source_DNA_AA_other_difference'])))
 chk(truth(d['at_least2_each_new_group_eligible'])==(truth(d['primary_eligible']) and int(d['primary_reused_guides'])>=2 and int(d['primary_fresh_guides'])>=2))
def aggregate(i,c,flag='primary_eligible',pc=1.):
 a=[oldv[(q['sgrna'],c)] for q in oldby[i] if truth(q[flag])]
 out={'old_mean':mean(a) if a else float('nan'),'old_median':statistics.median(a) if a else float('nan'),'old_guides':len(a)}
 for g in ('reused','fresh'):
  z=[q for q in newby[i] if truth(q[flag]) and q['guide_group']==g]
  vals=[float(newout[(q['sgrna'],c)]['LFC']) if pc==1 else newv[(q['sgrna'],c,pc)] for q in z]
  out[g+'_LFC']=mean(vals) if vals else float('nan');out[g+'_guides']=len(z);out[g+'_GC']=mean([float(q['GC_fraction_20nt']) for q in z]) if z else float('nan')
 out['fresh_minus_reused']=out['fresh_LFC']-out['reused_LFC']
 return out
arch={}
for fname,flag in [('all381_primary_product_cell_guide_group_statistics.csv','primary_eligible'),('all381_unique_only_product_cell_guide_group_statistics.csv','unique_only_eligible')]:
 vals=rows(O/fname);chk(len(vals)==381)
 for r in vals:
  v=aggregate(r['pair_ID'],r['cell'],flag)
  for k in v:
   if r[k]=='':chk(math.isnan(v[k]) or (k.endswith('_guides') and v[k]==0))
   else:eq(r[k],v[k])
 arch[flag]={(r['pair_ID'],r['cell']):r for r in vals}
pids=[r['pair_ID'] for r in den if truth(r['primary_eligible'])]
def arrays(ids,flag='primary_eligible',oldmetric='old_median',pc=1.):
 aa=[];rr=[];ff=[]
 for i in ids:
  av=[];rv=[];fv=[]
  for c in cells:
   a=arch[flag][(i,c)] if pc==1 else aggregate(i,c,flag,pc)
   av.append(float(a[oldmetric]));rv.append(float(a['reused_LFC']));fv.append(float(a['fresh_LFC']))
  aa.append(av);rr.append(rv);ff.append(fv)
 return tuple(np.array(x) for x in (aa,rr,ff))
A,R,F=arrays(pids);chk(len(pids)==88)
variants={}
for r in rows(O/'all_declared_primary_and_sensitivity_point_estimates.csv'):
 name=r['variant'];ids=pids;flag='primary_eligible';om='old_median';pc=1.
 if name=='older_mean':om='old_mean'
 if name=='unique_only_no_CDS_exclusion':flag='unique_only_eligible';ids=[d['pair_ID'] for d in den if truth(d[flag])]
 for nm,key in [('prior_source_AA_clean','source_AA_clean_eligible'),('expression_qualified_all3','expression_qualified_cohort_eligible'),('at_least2_each_new_group','at_least2_each_new_group_eligible')]:
  if name==nm:ids=[i for i in pids if truth(dm[i][key])]
 if name=='newer_pc0_5':pc=.5
 if name=='newer_pc5':pc=5.
 a,b,c=arrays(ids,flag,om,pc);v=stats(a,b,c);eq(r['fresh_minus_reused_LFC'],v[0]);eq(r['reused_minus_fresh_rank_agreement'],v[1]);chk(int(r['products'])==len(ids));variants[name]=(a,b,c)
for r in rows(O/'all_declared_cell_specific_cross_screen_rank_correlations.csv'):
 a,b,c=variants[r['variant']];k=cells.index(r['cell']);eq(r['rho_old_vs_reused'],rho(a[:,k],b[:,k]));eq(r['rho_old_vs_fresh'],rho(a[:,k],c[:,k]))
for r in rows(O/'all3_leave_one_common_cell_out_point_estimates.csv'):
 v=stats(A,R,F,tuple(i for i,c in enumerate(cells) if c!=r['omitted_cell']))
 eq(r['fresh_minus_reused_LFC'],v[0]);eq(r['reused_minus_fresh_rank_agreement'],v[1])
bp=json.loads(gzip.open(O/'all2000_primary_cluster_resampling_plans.json.gz','rt').read());bs=rows(O/'all2000_primary_cluster_bootstrap_statistics.csv');chk(bp['primary_pair_IDs']==pids)
rng=np.random.default_rng(20261015);expected=rng.integers(0,len(bp['clusters']),size=(2000,len(bp['clusters'])));chk(np.array_equal(expected,np.array(bp['draws'])))
for r,draw in zip(bs,bp['draws']):
 ix=[j for d in draw for j,i in enumerate(pids) if dm[i]['dependence_cluster']==bp['clusters'][d]]
 a,b,c=A[ix],R[ix],F[ix];v=stats(a,b,c);eq(r['fresh_minus_reused_LFC'],v[0]);eq(r['reused_minus_fresh_rank_agreement'],v[1])
 for k,cell in enumerate(cells):eq(r[cell+'_rho_old_vs_reused'],rho(a[:,k],b[:,k]));eq(r[cell+'_rho_old_vs_fresh'],rho(a[:,k],c[:,k]))
result=json.loads((O/'FROZEN_FRESH_GUIDE_VALIDATION_RESULTS.json').read_text())
def qlinear(v,p):
 v=sorted(v);x=(len(v)-1)*p;j=int(x);return v[j]+(v[min(j+1,len(v)-1)]-v[j])*(x-j)
for key,info in result['primary_intervals'].items():
 vals=[float(r[key]) for r in bs]
 for a,b in zip(info['CI97_5'],[qlinear(vals,.0125),qlinear(vals,.9875)]):eq(a,b)
 chk(info['undefined_draws']==0)
bplan=json.loads(gzip.open(O/'all1000_balanced_guide_selection_plans.json.gz','rt').read());bstat=rows(O/'all1000_balanced_guide_count_draw_statistics.csv')
vectors=[]
for i,inv in zip(pids,bplan['guide_inventory']):
 chk(inv['pair_ID']==i)
 groups=[]
 for g in ('reused','fresh'):
  guides=sorted(q['sgrna'] for q in newby[i] if truth(q['primary_eligible']) and q['guide_group']==g);chk(inv[g+'_guides']==guides)
  groups.append(np.array([[newv[(q,c,1.)] for c in cells] for q in guides]))
 chk(inv['balanced_count_each']==min(map(len,groups)));vectors.append((*groups,inv['balanced_count_each']))
for j,(choices,r) in enumerate(zip(bplan['draw_choices'],bstat)):
 rs=[];fs=[]
 for (br,bf,k),(ir,ifr) in zip(vectors,choices):
  chk(np.array_equal(rng.choice(len(br),size=k,replace=False),ir));chk(np.array_equal(rng.choice(len(bf),size=k,replace=False),ifr))
  rs.append(br[ir].mean(axis=0));fs.append(bf[ifr].mean(axis=0))
 b,c=np.array(rs),np.array(fs);v=stats(A,b,c);eq(r['fresh_minus_reused_LFC'],v[0]);eq(r['reused_minus_fresh_rank_agreement'],v[1])
 for k,cell in enumerate(cells):eq(r[cell+'_rho_old_vs_reused'],rho(A[:,k],b[:,k]));eq(r[cell+'_rho_old_vs_fresh'],rho(A[:,k],c[:,k]))
z=[arch['primary_eligible'][(i,cells[0])] for i in pids]
X=np.array([[1.,float(r['fresh_GC'])-float(r['reused_GC']),math.log2(float(r['fresh_guides'])/float(r['reused_guides']))] for r in z]);y=(F-R).mean(axis=1)
Q,T=np.linalg.qr(X,mode='reduced');coef=np.linalg.solve(T,Q.T@y)
adj=json.loads((O/'DESCRIPTIVE_GC_AND_GUIDECOUNT_ADJUSTED_WITHIN_ASSAY_DIFFERENCE.json').read_text())
for a,b in zip(coef,adj['coefficients']):eq(a,b)
report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass','checks':checks,'maximum_numeric_absolute_difference':maxerr,'source_guide_joins_roles_geometry_coverage_expression_and_dependence_checked':True,'all381_primary_and381_unique_only_product_cells_checked':True,'all8_point_variants_all24_cell_correlations_all3_omissions_checked':True,'all2000_primary_bootstrap_plans_and_statistics_checked':True,'all1000_balanced_guide_choices_and_statistics_checked':True,'manual_midrank_correlation_independent_of_scipy':True,'scalar_source_logarithm_checks_before_retaining_archived_float_ties':True,'QR_adjustment_independent_of_lstsq':True,'all10_source_SHA256_pins_checked':True}
(O/'INDEPENDENT_FRESH_GUIDE_SOURCE_AND_ALL_RESAMPLING_VERIFICATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
