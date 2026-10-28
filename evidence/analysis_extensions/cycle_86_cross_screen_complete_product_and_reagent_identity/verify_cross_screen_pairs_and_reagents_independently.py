from pathlib import Path
import csv,re,json,collections,datetime,os
O=Path(__file__).resolve().parent;N=O.parent;C66=N/'cycle_66_independent_nORF_fitness_benchmark';C34=N/'cycle_34_endogenous_peptide_and_ORF_provenance';C83=N/'cycle_83_human_nORF_guide_splice_liability'
def read(p):
 with p.open(newline='',encoding='utf8') as f:return list(csv.DictReader(f))
old=read(O/'all553_original2021_complete_sequence_models_before_cross_screen_comparison.csv');new=read(O/'all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv');assert len(old)==553 and len(new)==2196
oldmap={r['source_Excel_row']:r for r in old};newmap={r['source_ORF_ID']:r for r in new};oldname={r['source_name']:r for r in old}
for r in old+new:
 r['_mapped']=r['reference_mapped']=='True';r['_blocks']=[(int(a),int(b)) for a,b in re.findall(r'(\d+)-(\d+)',r['blocks_in_RNA_order0'])];r['_bases']={};r['_phase']={};n=0
 for j,(a,b) in enumerate(r['_blocks']):
  for base in range(a,b) if r['strand']=='+' else range(b-1,a-1,-1):r['_bases'][base]=j;r['_phase'][base]=n%3;n+=1
 assert len(r['_bases'])==int(float(r['DNA_nt'])) if r['_mapped'] else not r['_bases']
index=collections.defaultdict(list)
for r in new:
 for base,phase in r['_phase'].items():index[(r['chrom'],r['strand'],base)].append((r['source_ORF_ID'],phase))
overlap=collections.defaultdict(lambda:[0,0])
for r in old:
 for base,phase in r['_phase'].items():
  for newID,np in index.get((r['chrom'],r['strand'],base),[]):z=overlap[(r['source_Excel_row'],newID)];z[0]+=1;z[1]+=phase==np
def gene(s):return re.sub(r'\.\d+$','',s) if re.fullmatch(r'ENSG\d+(?:\.\d+)?',s) else ''
genes=collections.defaultdict(list)
for r in new:genes[gene(r['source_gene_ID_literal'])].append(r['source_ORF_ID'])
universe=set(overlap)
for r in old:
 g=gene(r['source_gene_ID_literal'])
 if g:universe.update((r['source_Excel_row'],v) for v in genes.get(g,[]))
prod=read(O/'all_source_gene_or_coding_overlap_cross_screen_product_pairs.csv');pm={(r['older_source_Excel_row'],r['newer_source_ORF_ID']):r for r in prod};assert len(pm)==len(prod)==274 and set(pm)==universe
# Reparse the original libraries and validated source-name crosswalk independently.
cross={r['source_S21_gene_name']:r['source_S2_name'] for r in read(C66/'all553_source_S2_S21_S23_gene_identifier_crosswalk.csv')}
olib={r['sgRNA Sequence']:cross[r['Gene_name']] for r in read(C66/'all_original_primary_CRISPR_guide_designs.csv') if r['Gene_name'] in cross}
nlib={r['sgrna']:r['gene'] for r in read(C34/'table13_guide_target_mapping.csv') if r['gene'] in newmap}
assert len(olib)==4391 and len(nlib)==12527
expectedreuse={(q,olib[q],nlib[q]) for q in olib.keys()&nlib.keys()};observedreuse={(r['sgrna'],r['older_source_name'],r['newer_source_ORF_ID']) for r in read(O/'all_exact20nt_guides_reused_across_source_ORF_libraries.csv')};assert expectedreuse==observedreuse and len(expectedreuse)==744
def inside(r,chrom,cut):return r['_mapped'] and chrom==r['chrom'] and cut-1 in r['_bases'] and cut in r['_bases'] and r['_bases'][cut-1]==r['_bases'][cut]
primary=re.compile(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)')
qual={};rawguides={}
for tag,folder,stem,lib,models in [('old',C83,'all_original2021_source_guide_exact',olib,oldname),('new',C34,'all_original_source_guide_exact',nlib,newmap)]:
 sites=collections.defaultdict(list)
 for r in read(folder/(stem+'_hg38_sites.csv')):
  q=r['sgrna']
  if q not in lib:continue
  assert r['source_target']==lib[q]
  pam=r['PAM_in_guide_orientation']
  if primary.fullmatch(r['chrom']) and len(pam)==3 and pam[1:]=='GG':sites[q].append((r['chrom'],int(r['cut_interbase0'])))
 qmodels=collections.defaultdict(dict);raw=collections.defaultdict(set)
 for q,target in lib.items():
  raw[target].add(q)
  if len(sites[q])==1:
   chrom,cut=sites[q][0]
   if inside(models[target],chrom,cut):qmodels[target][q]=(chrom,cut)
 qual[tag]=qmodels;rawguides[tag]=raw
assert sum(map(len,qual['old'].values()))==4169 and sum(map(len,qual['new'].values()))==12156
checks=0
def check(row,k,value):
 global checks
 got=row[k]
 if value is None:assert got=='',(k,got,value)
 elif isinstance(value,bool):assert got==str(value),(k,got,value)
 elif isinstance(value,int):assert float(got)==value,(k,got,value)
 else:assert got==str(value),(k,got,value)
 checks+=1
def stop(s):return s[:-1] if s.endswith('*') else s
classcounter=collections.Counter();expected={}
for key in sorted(universe):
 a,b=key;x=oldmap[a];y=newmap[b];r=pm[key];mapped=x['_mapped'] and y['_mapped'];nt,phase=overlap.get(key,[0,0]) if mapped else [None,None];samegene=bool(gene(x['source_gene_ID_literal']) and gene(x['source_gene_ID_literal'])==gene(y['source_gene_ID_literal']));chain=mapped and x['chrom']==y['chrom'] and x['strand']==y['strand'] and x['_blocks']==y['_blocks'];dna=mapped and x['DNA_sequence_in_source_stop_convention']==y['DNA_sequence_in_source_stop_convention'];aa=x['literal_standard_AA'];bb=y['literal_standard_AA'];full=mapped and aa==bb;start=mapped and x['chrom']==y['chrom'] and x['strand']==y['strand'] and float(x['coding_first_base0'])==float(y['coding_first_base0'])
 if not mapped:kind='related_gene_with_unresolved_older_reference'
 elif chain and dna:kind='exact_complete_reference_coding_chain_and_DNA'
 elif full:kind='exact_literal_complete_AA_with_different_chain_or_DNA'
 elif start:kind='same_genomic_start_nonidentical_complete_AA'
 elif phase:kind='different_start_with_same_phase_body_overlap'
 elif nt:kind='body_overlap_without_same_phase_coding_body'
 else:kind='same_source_gene_without_coding_body_overlap'
 classcounter[kind]+=1;og=rawguides['old'][x['source_name']];ng=rawguides['new'][b];oq=qual['old'][x['source_name']];nq=qual['new'][b];printed=stop(x['source_printed_AA_literal']);structure=stop(y['structure_full_AA_literal']);contains=(aa in bb or bb in aa) if mapped else None
 vals={'older_source_name':x['source_name'],'same_single_literal_stable_source_gene_ID':samegene,'older_reference_mapped':x['_mapped'],'same_strand_coding_body_overlap_nt':nt,'same_phase_body_overlap_nt':phase,'identity_relation_class':kind,'exact_complete_reference_coding_chain':chain if mapped else None,'exact_complete_mature_coding_DNA':dna if mapped else None,'exact_literal_complete_AA':full if mapped else None,'same_genomic_first_coding_base':start if mapped else None,'older_literal_AA_length':len(aa) if mapped else None,'newer_literal_AA_length':len(bb),'older_printed_source_AA_length':len(printed),'newer_structure_AA_length':len(structure),'exact_older_printed_AA_to_newer_structure_reference_AA':bool(printed==structure and structure),'older_source_printed_vs_literal_translation_agree':printed==aa if mapped else None,'newer_structure_vs_literal_translation_agree':structure==bb,'exact_initiator_M_sensitivity_full_AA':x['initiator_M_sensitivity_AA']==y['initiator_M_sensitivity_AA'] if mapped else None,'complete_AA_one_contains_other':contains,'AA_common_prefix_length':len(os.path.commonprefix([aa,bb])) if mapped else None,'AA_common_suffix_length':len(os.path.commonprefix([aa[::-1],bb[::-1]])) if mapped else None,'older_assigned_guide_count':len(og),'newer_assigned_guide_count':len(ng),'exact20nt_reused_guide_count':len(og&ng),'older_qualifying_unique_primary_NGG_inside_own_ORF_guides':len(oq),'newer_qualifying_unique_primary_NGG_inside_own_ORF_guides':len(nq),'reused_guides_qualifying_in_both_ORFs':len(oq.keys()&nq.keys()),'older_qualifying_guides_cut_inside_both_bodies':sum(inside(y,ch,c) for ch,c in oq.values()) if mapped else None,'older_qualifying_guides_cut_only_older_body':sum(not inside(y,ch,c) for ch,c in oq.values()) if mapped else None,'newer_qualifying_guides_cut_inside_both_bodies':sum(inside(x,ch,c) for ch,c in nq.values()) if mapped else None,'newer_qualifying_guides_cut_only_newer_body':sum(not inside(x,ch,c) for ch,c in nq.values()) if mapped else None}
 for k,v in vals.items():check(r,k,v)
 expected[key]=vals
for r in read(O/'all553_original2021_cross_screen_identity_denominators.csv'):
 a=r['source_Excel_row'];z=[v for (i,j),v in expected.items() if i==a]
 for k,v in {'related_native_targets':len(z),'same_source_gene_native_targets':sum(v['same_single_literal_stable_source_gene_ID'] for v in z),'full_reference_chain_matches':sum(v['exact_complete_reference_coding_chain'] is True for v in z),'literal_complete_AA_matches':sum(v['exact_literal_complete_AA'] is True for v in z),'printed_source_AA_to_structure_exact_matches':sum(v['exact_older_printed_AA_to_newer_structure_reference_AA'] for v in z),'same_start_nonidentical_product_targets':sum(v['identity_relation_class']=='same_genomic_start_nonidentical_complete_AA' for v in z),'related_targets_reusing_at_least_one_guide':sum(v['exact20nt_reused_guide_count']>0 for v in z)}.items():check(r,k,v)
for r in read(O/'all2196_original2026_cross_screen_identity_denominators.csv'):
 b=r['source_ORF_ID'];z=[v for (i,j),v in expected.items() if j==b]
 for k,v in {'related_older_targets':len(z),'same_source_gene_older_targets':sum(v['same_single_literal_stable_source_gene_ID'] for v in z),'full_reference_chain_matches':sum(v['exact_complete_reference_coding_chain'] is True for v in z),'literal_complete_AA_matches':sum(v['exact_literal_complete_AA'] is True for v in z),'printed_source_AA_to_structure_exact_matches':sum(v['exact_older_printed_AA_to_newer_structure_reference_AA'] for v in z),'same_start_nonidentical_product_targets':sum(v['identity_relation_class']=='same_genomic_start_nonidentical_complete_AA' for v in z),'related_targets_reusing_at_least_one_guide':sum(v['exact20nt_reused_guide_count']>0 for v in z)}.items():check(r,k,v)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass','checks':checks,'all274_related_product_pairs_reconstructed_by_per_base_frame_index':True,'all553_and2196_denominator_records_checked':True,'all4391_and12527_original_guide_assignments_reparsed':True,'all744_exact_reused_guide_strings_verified':True,'qualified_guides_reconstructed_from_primary_chromosome_and_literal_PAM_not_stored_flags':{'older':4169,'newer':12156},'all_pair_cut_membership_fields_checked_with_chromosome_and_exon_membership':True,'identity_class_counts':dict(classcounter),'no_functional_outcome_join_or_new_genome_scan':True};(O/'INDEPENDENT_CROSS_SCREEN_PAIR_AND_REAGENT_VERIFICATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
