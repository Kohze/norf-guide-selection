from pathlib import Path
import json,itertools,hashlib,datetime
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;N=O.parent
C70=N/'cycle_70_ORF_resolved_context_and_host_dependency';C86=N/'cycle_86_cross_screen_complete_product_and_reagent_identity';C83=N/'cycle_83_human_nORF_guide_splice_liability';C35=N/'cycle_35_guide_replicates_and_medical_provenance'
read=lambda p,**kw:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',low_memory=False,**kw)
meta=read(C70/'all17568_source_resolved_context_host_fitness_and_eligibility.csv',usecols=['cell_line','orf_targeted','Ensembl_id','riboseq_tpm','rnaseq_tpm','dependence_union_cluster'])
models=read(C86/'all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv').set_index('source_ORF_ID')
geom=read(C83/'all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv')
counts=read(C35/'all_original_table13_B_guide_measurements.csv',usecols=['sgrna','pDNA_1','pDNA_2','pDNA_3']).set_index('sgrna')
oldmatch=read(C86/'all_source_gene_or_coding_overlap_cross_screen_product_pairs.csv');oldmatch=oldmatch[oldmatch.identity_relation_class.eq('exact_complete_reference_coding_chain_and_DNA')]
oldmap=dict(zip(oldmatch.newer_source_ORF_ID,oldmatch.older_source_name))
truth=lambda x:str(x).lower()=='true'
native=meta.drop_duplicates('orf_targeted');assert len(native)==2196
guideby={k:v.to_dict('records') for k,v in geom.groupby('source_target')}
rows=[];guide_rows=[];pairs=[]
for gene,g in native[native.Ensembl_id.str.match(r'^ENSG[0-9]+$')].groupby('Ensembl_id',sort=True):
 ids=sorted(g.orf_targeted)
 for a,b in itertools.combinations(ids,2):
  pair=gene+'__'+a+'__'+b;pairs.append(pair);out={'pair_ID':pair,'source_gene':gene,'ORF_A':a,'ORF_B':b,'biotype_A':models.loc[a,'source_biotype'],'biotype_B':models.loc[b,'source_biotype'],'both_exact_old_panel_products':a in oldmap and b in oldmap,'old_source_A':oldmap.get(a,''),'old_source_B':oldmap.get(b,'')}
  clusters=set(meta.loc[meta.orf_targeted.isin([a,b]),'dependence_union_cluster']);assert len(clusters)==1;out['dependence_cluster']=next(iter(clusters))
  for role,own,partner in [('A',a,b),('B',b,a)]:
   m=models.loc[partner];blocks=[tuple(map(int,s.split('-'))) for s in m.blocks_in_RNA_order0.split('|')];good=[];safe=[]
   for q in guideby.get(own,[]):
    base=truth(q['primary_geometry_eligible']);coverage=bool((counts.loc[q['sgrna']].to_numpy(float)>=10).all())
    intersects=False
    if base and q['chrom']==m.chrom:
     lo=int(float(q['protospacer_start0']));hi=int(float(q['protospacer_end0']));intersects=any(lo<int(y) and hi>int(x) for x,y in blocks)
    eligible=base and coverage and not intersects;is_safe=eligible and truth(q['splice_safe_gt50'])
    guide_rows.append({'pair_ID':pair,'source_gene':gene,'role':role,'own_ORF':own,'partner_ORF':partner,'sgrna':q['sgrna'],'primary_geometry_eligible':base,'plasmid_coverage_pass':coverage,'protospacer_intersects_partner_coding_body':intersects,'pair_specific_eligible':eligible,'CDS_bearing_RNA_splice_safe_gt50':is_safe,'GC_fraction':q['GC_fraction_20nt']})
    if eligible:good.append(q['sgrna'])
    if is_safe:safe.append(q['sgrna'])
   out[role+'_guide_count']=len(good);out[role+'_guides']=json.dumps(sorted(good));out[role+'_splice_safe_count']=len(safe);out[role+'_splice_safe_guides']=json.dumps(sorted(safe))
  z=meta[meta.orf_targeted.isin([a,b])];allfinite=True
  for v in ['riboseq_tpm','rnaseq_tpm']:
   x=pd.to_numeric(z[v],errors='coerce');allfinite=allfinite and np.isfinite(x).all()
  assert len(z)==16
  out['all8_expression_context_finite']=bool(allfinite);out['primary_metadata_eligible']=allfinite and out['A_guide_count']>=4 and out['B_guide_count']>=4;out['min3_metadata_eligible']=allfinite and out['A_guide_count']>=3 and out['B_guide_count']>=3;out['splice_safe_metadata_eligible']=allfinite and out['A_splice_safe_count']>=4 and out['B_splice_safe_count']>=4
  rows.append(out)
den=pd.DataFrame(rows);assert len(den)==145 and den.source_gene.nunique()==117
den.to_csv(O/'all145_within_gene_pair_eligibility_before_outcome_join.csv',index=False)
pd.DataFrame(guide_rows).to_csv(O/'all_within_gene_pair_specific_guide_eligibility_before_outcome_join.csv',index=False)
meta[meta.orf_targeted.isin(set(den.ORF_A)|set(den.ORF_B))].to_csv(O/'all1984_multiORF_expression_context_before_outcome_join.csv',index=False)
inputs=[C70/'all17568_source_resolved_context_host_fitness_and_eligibility.csv',C86/'all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv',C83/'all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv',C35/'all_original_table13_B_guide_measurements.csv',C86/'all_source_gene_or_coding_overlap_cross_screen_product_pairs.csv',C83/'all215704_original2026_guide_cell_primary_log_abundance_changes.csv']
pins=[{'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs]
(O/'PRE_OUTCOME_INPUT_PINS.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pins':pins},indent=2))
cohorts={}
for flag in ['primary_metadata_eligible','min3_metadata_eligible','splice_safe_metadata_eligible']:
 z=den[den[flag]];cohorts[flag]={'pairs':len(z),'genes':z.source_gene.nunique(),'clusters':z.dependence_cluster.nunique(),'both_exact_old_product_pairs':int(z.both_exact_old_panel_products.sum())}
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'metadata_preflight_complete_before_outcomes','all_pairs':145,'all_source_genes':117,'cohorts':cohorts,'all_pairs_both_exact_old_products':int(den.both_exact_old_panel_products.sum()),'guide_assignment_rows':len(guide_rows),'no_new_pair_conditioned_outcomes_inspected':True}
(O/'WITHIN_LOCUS_METADATA_FEASIBILITY_RESULTS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
