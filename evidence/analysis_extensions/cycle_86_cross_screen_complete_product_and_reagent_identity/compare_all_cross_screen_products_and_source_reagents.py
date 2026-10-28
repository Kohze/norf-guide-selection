from pathlib import Path
import pandas as pd,numpy as np,json,re,datetime,collections
O=Path(__file__).resolve().parent;N=O.parent;C34=N/'cycle_34_endogenous_peptide_and_ORF_provenance';C83=N/'cycle_83_human_nORF_guide_splice_liability'
old=pd.read_csv(O/'all553_original2021_complete_sequence_models_before_cross_screen_comparison.csv',keep_default_na=False);new=pd.read_csv(O/'all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv',keep_default_na=False);oldrows=old.to_dict('records');newrows=new.to_dict('records')
truth=lambda v:str(v).lower() in ['true','1'];strip=lambda s:s[:-1] if s.endswith('*') else s
for r in oldrows+newrows:
 r['mapped']=truth(r['reference_mapped']);
 if r['mapped']:
  coordinate=float(r['coding_first_base0']);assert coordinate.is_integer();r['coding_first_base0']=int(coordinate)
 r['blocks']=[(int(a),int(b)) for a,b in re.findall(r'(\d+)-(\d+)',r.get('blocks_in_RNA_order0',''))];r['offsets']=[];offset=0
 for a,b in r['blocks']:r['offsets'].append(offset);offset+=b-a
 r['body_bounding_lo']=min((a for a,b in r['blocks']),default=None);r['body_bounding_hi']=max((b for a,b in r['blocks']),default=None)
newgene=collections.defaultdict(list);newloc=collections.defaultdict(list)
for i,r in enumerate(newrows):
 if r['source_single_stable_gene_ID']:newgene[r['source_single_stable_gene_ID']].append(i)
 newloc[(r['chrom'],r['strand'])].append(i)
def overlap(x,y):
 total=phase=0
 if not x['mapped'] or not y['mapped'] or x['chrom']!=y['chrom'] or x['strand']!=y['strand']:return None,None
 for i,(a,b) in enumerate(x['blocks']):
  for j,(c,d) in enumerate(y['blocks']):
   lo,hi=max(a,c),min(b,d)
   if lo>=hi:continue
   total+=hi-lo
   p=(x['offsets'][i]+(lo-a if x['strand']=='+' else b-1-lo))%3;q=(y['offsets'][j]+(lo-c if y['strand']=='+' else d-1-lo))%3
   if p==q:phase+=hi-lo
 return total,phase
def shared_prefix(x,y):return next((i for i,(a,b) in enumerate(zip(x,y)) if a!=b),min(len(x),len(y)))
def inside(r,cut):return any(a<cut<b for a,b in r['blocks'])
guides={};qualified={};guide_pins=[]
for cohort,folder,names in [('older',C83,{r['source_name'] for r in oldrows}),('newer',C34,{r['source_ORF_ID'] for r in newrows})]:
 stem='all_original2021_source_guide_exact' if cohort=='older' else 'all_original_source_guide_exact';summ=pd.read_csv(folder/(stem+'_genome_summary.csv'),keep_default_na=False);summ=summ[summ.source_target.isin(names)];assert not summ.sgrna.duplicated().any();guides[cohort]={name:set(g.sgrna) for name,g in summ.groupby('source_target')};good=summ[(summ.primary_NGG_sites.astype(int)==1)&(summ.primary_NGG_cut_in_assigned_ORF.astype(int)==1)];sites=pd.read_csv(folder/(stem+'_hg38_sites.csv'),usecols=['sgrna','source_target','chrom','primary_chromosome','canonical_NGG_PAM','cut_interbase0'],keep_default_na=False);sites=sites[sites.primary_chromosome.map(truth)&sites.canonical_NGG_PAM.map(truth)];z=sites.merge(good[['sgrna','source_target']],on=['sgrna','source_target'],how='inner',validate='many_to_one');assert len(z)==len(good) and not z[['sgrna','source_target']].duplicated().any()
 qualified[cohort]={name:dict(zip(g.sgrna,g.cut_interbase0.astype(int))) for name,g in z.groupby('source_target')}
 if cohort=='older':oldguide=summ[['sgrna','source_target']].rename(columns={'source_target':'older_source_name'})
 else:newguide=summ[['sgrna','source_target']].rename(columns={'source_target':'newer_source_ORF_ID'})
reuse=oldguide.merge(newguide,on='sgrna',how='inner',validate='one_to_one');reuse.to_csv(O/'all_exact20nt_guides_reused_across_source_ORF_libraries.csv',index=False)
pairs=[]
for x in oldrows:
 candidates=set(newgene.get(x['source_single_stable_gene_ID'],[])) if x['source_single_stable_gene_ID'] else set()
 if x['mapped']:
  for j in newloc.get((x['chrom'],x['strand']),[]):
   y=newrows[j]
   if x['body_bounding_lo']<y['body_bounding_hi'] and y['body_bounding_lo']<x['body_bounding_hi']:candidates.add(j)
 for j in sorted(candidates):
  y=newrows[j];gene=bool(x['source_single_stable_gene_ID'] and x['source_single_stable_gene_ID']==y['source_single_stable_gene_ID']);nt,phase=overlap(x,y)
  if not gene and not(nt and nt>0):continue
  xAA=x.get('literal_standard_AA','');yAA=y.get('literal_standard_AA','');mapped=x['mapped'] and y['mapped'];chain=bool(mapped and x['chrom']==y['chrom'] and x['strand']==y['strand'] and x['blocks']==y['blocks']);dna=bool(mapped and x['DNA_sequence_in_source_stop_convention']==y['DNA_sequence_in_source_stop_convention']);fullAA=bool(mapped and xAA==yAA);start=bool(mapped and x['chrom']==y['chrom'] and x['strand']==y['strand'] and int(x['coding_first_base0'])==int(y['coding_first_base0']));structureAA=strip(y.get('structure_full_AA_literal',''));printedAA=x['source_printed_AA_without_one_terminal_stop'];sourceStructure=bool(structureAA and printedAA==structureAA)
  if not mapped:kind='related_gene_with_unresolved_older_reference'
  elif chain and dna:kind='exact_complete_reference_coding_chain_and_DNA'
  elif fullAA:kind='exact_literal_complete_AA_with_different_chain_or_DNA'
  elif start:kind='same_genomic_start_nonidentical_complete_AA'
  elif phase:kind='different_start_with_same_phase_body_overlap'
  elif nt:kind='body_overlap_without_same_phase_coding_body'
  else:kind='same_source_gene_without_coding_body_overlap'
  og=guides['older'].get(x['source_name'],set());ng=guides['newer'].get(y['source_ORF_ID'],set());oq=qualified['older'].get(x['source_name'],{});nq=qualified['newer'].get(y['source_ORF_ID'],{});reused=og&ng
  row={'older_source_Excel_row':x['source_Excel_row'],'older_source_name':x['source_name'],'newer_source_Excel_row':y['source_Excel_row'],'newer_source_ORF_ID':y['source_ORF_ID'],'older_source_gene_ID':x['source_gene_ID_literal'],'newer_source_gene_ID':y['source_gene_ID_literal'],'same_single_literal_stable_source_gene_ID':gene,'older_reference_mapped':x['mapped'],'same_strand_coding_body_overlap_nt':nt,'same_phase_body_overlap_nt':phase,'identity_relation_class':kind,'exact_complete_reference_coding_chain':chain if mapped else None,'exact_complete_mature_coding_DNA':dna if mapped else None,'exact_literal_complete_AA':fullAA if mapped else None,'same_genomic_first_coding_base':start if mapped else None,'older_literal_AA_length':len(xAA) if mapped else None,'newer_literal_AA_length':len(yAA),'older_printed_source_AA_length':len(printedAA),'newer_structure_AA_length':len(structureAA),'exact_older_printed_AA_to_newer_structure_reference_AA':sourceStructure,'older_source_printed_vs_literal_translation_agree':printedAA==xAA if mapped else None,'newer_structure_vs_literal_translation_agree':structureAA==yAA,'older_prior_source_AA_other_difference':x['source_AA_other_difference_prior_flag'],'exact_initiator_M_sensitivity_full_AA':x.get('initiator_M_sensitivity_AA')==y.get('initiator_M_sensitivity_AA') if mapped else None,'complete_AA_one_contains_other':bool(xAA in yAA or yAA in xAA) if mapped else None,'AA_common_prefix_length':shared_prefix(xAA,yAA) if mapped else None,'AA_common_suffix_length':shared_prefix(xAA[::-1],yAA[::-1]) if mapped else None,'older_assigned_guide_count':len(og),'newer_assigned_guide_count':len(ng),'exact20nt_reused_guide_count':len(reused),'older_qualifying_unique_primary_NGG_inside_own_ORF_guides':len(oq),'newer_qualifying_unique_primary_NGG_inside_own_ORF_guides':len(nq),'reused_guides_qualifying_in_both_ORFs':len(set(oq)&set(nq)),'older_qualifying_guides_cut_inside_both_bodies':sum(x['chrom']==y['chrom'] and inside(y,cut) for cut in oq.values()) if mapped else None,'older_qualifying_guides_cut_only_older_body':sum(x['chrom']!=y['chrom'] or not inside(y,cut) for cut in oq.values()) if mapped else None,'newer_qualifying_guides_cut_inside_both_bodies':sum(x['chrom']==y['chrom'] and inside(x,cut) for cut in nq.values()) if mapped else None,'newer_qualifying_guides_cut_only_newer_body':sum(x['chrom']!=y['chrom'] or not inside(x,cut) for cut in nq.values()) if mapped else None}
  for q,cut in oq.items():assert inside(x,cut)
  for q,cut in nq.items():assert inside(y,cut)
  pairs.append(row)
p=pd.DataFrame(pairs);assert not p[['older_source_Excel_row','newer_source_ORF_ID']].duplicated().any();p.to_csv(O/'all_source_gene_or_coding_overlap_cross_screen_product_pairs.csv',index=False)
old_den=[]
for x in oldrows:
 z=p[p.older_source_Excel_row.eq(x['source_Excel_row'])];old_den.append({'source_Excel_row':x['source_Excel_row'],'source_name':x['source_name'],'reference_mapped':x['mapped'],'related_native_targets':len(z),'same_source_gene_native_targets':int(z.same_single_literal_stable_source_gene_ID.sum()),'full_reference_chain_matches':int(z.exact_complete_reference_coding_chain.eq(True).sum()),'literal_complete_AA_matches':int(z.exact_literal_complete_AA.eq(True).sum()),'printed_source_AA_to_structure_exact_matches':int(z.exact_older_printed_AA_to_newer_structure_reference_AA.sum()),'same_start_nonidentical_product_targets':int(z.identity_relation_class.eq('same_genomic_start_nonidentical_complete_AA').sum()),'related_targets_reusing_at_least_one_guide':int(z.exact20nt_reused_guide_count.gt(0).sum())})
pd.DataFrame(old_den).to_csv(O/'all553_original2021_cross_screen_identity_denominators.csv',index=False)
new_den=[]
for y in newrows:
 z=p[p.newer_source_ORF_ID.eq(y['source_ORF_ID'])];new_den.append({'source_ORF_ID':y['source_ORF_ID'],'related_older_targets':len(z),'same_source_gene_older_targets':int(z.same_single_literal_stable_source_gene_ID.sum()),'full_reference_chain_matches':int(z.exact_complete_reference_coding_chain.eq(True).sum()),'literal_complete_AA_matches':int(z.exact_literal_complete_AA.eq(True).sum()),'printed_source_AA_to_structure_exact_matches':int(z.exact_older_printed_AA_to_newer_structure_reference_AA.sum()),'same_start_nonidentical_product_targets':int(z.identity_relation_class.eq('same_genomic_start_nonidentical_complete_AA').sum()),'related_targets_reusing_at_least_one_guide':int(z.exact20nt_reused_guide_count.gt(0).sum())})
pd.DataFrame(new_den).to_csv(O/'all2196_original2026_cross_screen_identity_denominators.csv',index=False)
classes=p.groupby('identity_relation_class').agg(pairs=('older_source_Excel_row','size'),older_ORFs=('older_source_Excel_row','nunique'),newer_ORFs=('newer_source_ORF_ID','nunique'),pairs_with_reused_guide=('exact20nt_reused_guide_count',lambda z:int((z>0).sum())),reused_pair_guide_assignments=('exact20nt_reused_guide_count','sum')).reset_index();classes.to_csv(O/'all_cross_screen_identity_class_and_guide_reuse_denominators.csv',index=False)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'complete_descriptive_identity_and_reagent_analysis_pending_independent_verification','all_older_candidates':553,'all_newer_targets':2196,'related_source_pairs':len(p),'older_ORFs_with_related_newer_targets':int(p.older_source_Excel_row.nunique()),'newer_ORFs_with_related_older_targets':int(p.newer_source_ORF_ID.nunique()),'identity_classes':classes.to_dict('records'),'older_assigned_unique_ORF_guides':len(oldguide),'newer_assigned_unique_ORF_guides':len(newguide),'exact20nt_reused_guide_strings':len(reuse),'older_guides_qualifying_unique_primary_NGG_inside_own_ORF':sum(len(v) for v in qualified['older'].values()),'newer_guides_qualifying_unique_primary_NGG_inside_own_ORF':sum(len(v) for v in qualified['newer'].values()),'older_ORFs_with_any_exact_full_reference_chain_match':sum(x['full_reference_chain_matches']>0 for x in old_den),'older_ORFs_with_any_exact_literal_full_AA_match':sum(x['literal_complete_AA_matches']>0 for x in old_den),'older_ORFs_with_any_source_printed_full_AA_to_structure_match':sum(x['printed_source_AA_to_structure_exact_matches']>0 for x in old_den),'no_functional_outcome_join_or_new_score_model':True,'medical_target_or_new_translation_claim':False};(O/'COMPLETE_CROSS_SCREEN_PRODUCT_AND_REAGENT_IDENTITY_RESULTS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
