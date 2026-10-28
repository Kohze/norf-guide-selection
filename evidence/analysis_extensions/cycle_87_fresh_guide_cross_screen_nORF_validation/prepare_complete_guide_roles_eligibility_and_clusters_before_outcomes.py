from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime
O=Path(__file__).resolve().parent;N=O.parent;C86=N/'cycle_86_cross_screen_complete_product_and_reagent_identity';C83=N/'cycle_83_human_nORF_guide_splice_liability';C66=N/'cycle_66_independent_nORF_fitness_benchmark';C70=N/'cycle_70_ORF_resolved_context_and_host_dependency';C34=N/'cycle_34_endogenous_peptide_and_ORF_provenance';C35=N/'cycle_35_guide_replicates_and_medical_provenance'
read=lambda p,**kw:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',**kw)
p=read(C86/'all_source_gene_or_coding_overlap_cross_screen_product_pairs.csv');p=p[p.identity_relation_class.eq('exact_complete_reference_coding_chain_and_DNA')].sort_values('older_source_Excel_row').copy();assert len(p)==127 and p.older_source_Excel_row.is_unique and p.newer_source_ORF_ID.is_unique;p['pair_ID']=['old'+str(x)+'__'+y for x,y in zip(p.older_source_Excel_row,p.newer_source_ORF_ID)]
oldlib=read(C83/'original2021_complete_guide_library_for_genome_scan.csv');assert len(oldlib)==5776 and oldlib.sgrna.is_unique;oldall=set(oldlib.sgrna);oldassigned={name:set(x.sgrna) for name,x in oldlib.groupby('gene')}
oldgeom=read(C83/'all_original2021_library_guide_splice_and_CDS_geometry_before_outcome_join.csv');newgeom=read(C83/'all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv');truth=lambda s:s.astype(str).str.lower().eq('true')
oldgeom['eligible_unique_primary_inside']=truth(oldgeom.unique_primary_NGG_site)&truth(oldgeom.cut_inside_assigned_source_ORF)&oldgeom.source_distinct_target_assignments_for_guide.eq(1);newgeom['eligible_unique_primary_inside']=truth(newgeom.unique_primary_NGG_site)&truth(newgeom.cut_inside_assigned_source_ORF)&newgeom.source_distinct_target_assignments_for_guide.eq(1)
presenceold=read(C66/'all_original_primary_screen_guide_LFC_rows.csv',usecols=['Cell Line','Construct Barcode']);presencenew=read(C83/'all215704_original2026_guide_cell_primary_log_abundance_changes.csv',usecols=['sgrna','cell','plasmid_coverage_pass']);cells=['A375','A549','HEPG2'];assert sorted(set(presenceold['Cell Line'])&set(presencenew.cell))==cells
coverage=presencenew.groupby('sgrna').plasmid_coverage_pass.nunique();assert coverage.eq(1).all();cover=presencenew.drop_duplicates('sgrna').set_index('sgrna').plasmid_coverage_pass.astype(str).str.lower().eq('true')
oldpresence=set(zip(presenceold['Construct Barcode'],presenceold['Cell Line']));newpresence=set(zip(presencenew.sgrna,presencenew.cell))
oldrows=[];newrows=[]
for r in p.to_dict('records'):
 a=oldgeom[oldgeom.source_target.eq(r['older_source_name'])];b=newgeom[newgeom.source_target.eq(r['newer_source_ORF_ID'])]
 assert len(a)==r['older_assigned_guide_count'] and len(b)==r['newer_assigned_guide_count']
 for q in a.to_dict('records'):
  oldrows.append({**q,'pair_ID':r['pair_ID'],'old_source_name':r['older_source_name'],'new_source_ORF_ID':r['newer_source_ORF_ID'],'source_present_all3sharedcells':all((q['sgrna'],c) in oldpresence for c in cells)})
 for q in b.to_dict('records'):
  group='reused' if q['sgrna'] in oldassigned[r['older_source_name']] else ('other_old_library_reuse' if q['sgrna'] in oldall else 'fresh')
  newrows.append({**q,'pair_ID':r['pair_ID'],'old_source_name':r['older_source_name'],'new_source_ORF_ID':r['newer_source_ORF_ID'],'guide_group':group,'plasmid_coverage_pass':bool(cover.loc[q['sgrna']]),'source_present_all3sharedcells':all((q['sgrna'],c) in newpresence for c in cells)})
og=pd.DataFrame(oldrows);ng=pd.DataFrame(newrows);assert len(og)==1010 and len(ng)==720
og['primary_eligible']=truth(og.primary_geometry_eligible)&og.source_present_all3sharedcells;og['unique_only_eligible']=og.eligible_unique_primary_inside&og.source_present_all3sharedcells
ng['primary_eligible']=truth(ng.primary_geometry_eligible)&ng.plasmid_coverage_pass&ng.source_present_all3sharedcells;ng['unique_only_eligible']=ng.eligible_unique_primary_inside&ng.plasmid_coverage_pass&ng.source_present_all3sharedcells
og.to_csv(O/'all1010_older_exact_product_guide_eligibility_before_outcome_join.csv',index=False);ng.to_csv(O/'all720_newer_exact_product_guide_roles_and_eligibility_before_outcome_join.csv',index=False)
expr=read(C34/'table13_new_eight_cell_fitness.csv',usecols=['orf_targeted','cell_line','riboseq_tpm','rnaseq_tpm']);expr=expr[expr.orf_targeted.isin(p.newer_source_ORF_ID)&expr.cell_line.isin(cells)].copy();assert len(expr)==127*3
expr['RNA_numeric']=pd.to_numeric(expr.rnaseq_tpm,errors='coerce');expr['Ribo_numeric']=pd.to_numeric(expr.riboseq_tpm,errors='coerce');expr['source_expression_known']=expr.RNA_numeric.notna()&expr.Ribo_numeric.notna();expr['source_expression_qualified']=expr.source_expression_known&expr.RNA_numeric.ge(10)&expr.Ribo_numeric.ge(5);expr.to_csv(O/'all381_source_expression_context_rows_before_outcome_join.csv',index=False)
oldclusters=read(C66/'all553_analysis_source_labels_context_and_locus_groups.csv',usecols=['source_Excel_row','source_locus_cluster']).set_index('source_Excel_row').source_locus_cluster.to_dict();newclusters=read(C70/'all17568_source_resolved_context_host_fitness_and_eligibility.csv',usecols=['orf_targeted','dependence_union_cluster']);assert newclusters.groupby('orf_targeted').dependence_union_cluster.nunique().eq(1).all();nc=newclusters.drop_duplicates('orf_targeted').set_index('orf_targeted').dependence_union_cluster.to_dict()
parent=list(range(len(p)));owners={}
def find(x):
 while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
 return x
def union(a,b):a,b=find(a),find(b);parent[max(a,b)]=min(a,b)
for i,r in enumerate(p.to_dict('records')):
 for token in [('older_locus',oldclusters[r['older_source_Excel_row']]),('newer_gene_locus',nc[r['newer_source_ORF_ID']])]:
  assert token[1]!=''
  if token in owners:union(i,owners[token])
  else:owners[token]=i
p['dependence_cluster']=['cluster'+str(find(i)).zfill(3) for i in range(len(p))]
rows=[]
for r in p.to_dict('records'):
 a=og[og.pair_ID.eq(r['pair_ID'])];b=ng[ng.pair_ID.eq(r['pair_ID'])];e=expr[expr.orf_targeted.eq(r['newer_source_ORF_ID'])];out={'pair_ID':r['pair_ID'],'older_source_Excel_row':r['older_source_Excel_row'],'older_source_name':r['older_source_name'],'newer_source_ORF_ID':r['newer_source_ORF_ID'],'dependence_cluster':r['dependence_cluster'],'prior_source_DNA_AA_other_difference':str(r['older_prior_source_AA_other_difference']).lower()=='true','expression_known_all3cells':bool(e.source_expression_known.all()),'expression_qualified_all3cells':bool(e.source_expression_qualified.all())}
 for tag,maskold,masknew in [('primary',a.primary_eligible,b.primary_eligible),('unique_only',a.unique_only_eligible,b.unique_only_eligible)]:
  oldn=int(maskold.sum());reuse=int((masknew&b.guide_group.eq('reused')).sum());fresh=int((masknew&b.guide_group.eq('fresh')).sum());out[tag+'_older_guides']=oldn;out[tag+'_reused_guides']=reuse;out[tag+'_fresh_guides']=fresh;out[tag+'_eligible']=oldn>=1 and reuse>=1 and fresh>=1
 out['at_least2_each_new_group_eligible']=out['primary_eligible'] and out['primary_reused_guides']>=2 and out['primary_fresh_guides']>=2;out['expression_qualified_cohort_eligible']=out['primary_eligible'] and out['expression_qualified_all3cells'];out['source_AA_clean_eligible']=out['primary_eligible'] and not out['prior_source_DNA_AA_other_difference'];rows.append(out)
den=pd.DataFrame(rows);den.to_csv(O/'all127_exact_product_cohort_eligibility_and_clusters_before_outcome_join.csv',index=False)
inputs=[C86/'all_source_gene_or_coding_overlap_cross_screen_product_pairs.csv',C83/'original2021_complete_guide_library_for_genome_scan.csv',C83/'all_original2021_library_guide_splice_and_CDS_geometry_before_outcome_join.csv',C83/'all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv',C66/'all_original_primary_screen_guide_LFC_rows.csv',C83/'all215704_original2026_guide_cell_primary_log_abundance_changes.csv',C35/'all_original_table13_B_guide_measurements.csv',C34/'table13_new_eight_cell_fitness.csv',C66/'all553_analysis_source_labels_context_and_locus_groups.csv',C70/'all17568_source_resolved_context_host_fitness_and_eligibility.csv']
pins=[{'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in inputs];(O/'PRE_OUTCOME_JOIN_INPUT_AND_SOURCE_PINS.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pins':pins},indent=2))
cohorts={k:{'products':int(den[k].sum()),'clusters':int(den.loc[den[k],'dependence_cluster'].nunique())} for k in ['primary_eligible','unique_only_eligible','at_least2_each_new_group_eligible','expression_qualified_cohort_eligible','source_AA_clean_eligible']}
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'complete_metadata_eligibility_and_dependence_preflight_before_new_outcome_joins','all_exact_product_pairs':127,'shared_cells':cells,'older_assigned_guides':len(og),'newer_assigned_guides':len(ng),'newer_guide_roles_all':ng.guide_group.value_counts().to_dict(),'newer_guide_roles_primary_eligible':ng.loc[ng.primary_eligible,'guide_group'].value_counts().to_dict(),'cohorts':cohorts,'primary_inferential_cluster_gate_satisfied':cohorts['primary_eligible']['clusters']>=10,'no_new_guide_group_effects_correlations_or_outcome_values_inspected':True};(O/'FRESH_GUIDE_ELIGIBILITY_AND_DEPENDENCE_PREFLIGHT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
