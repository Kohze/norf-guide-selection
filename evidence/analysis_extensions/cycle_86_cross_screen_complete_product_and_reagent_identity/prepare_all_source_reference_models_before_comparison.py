from pathlib import Path
import pandas as pd,numpy as np,json,re,csv,struct,itertools,hashlib,datetime,collections
O=Path(__file__).resolve().parent;N=O.parent;R=N.parent.parent;C66=N/'cycle_66_independent_nORF_fitness_benchmark';C34=N/'cycle_34_endogenous_peptide_and_ORF_provenance';REF=R/'full_database_study_2026-09-13/references/hg38.2bit'
sourcefiles=[C66/'all553_original_ORFeome_identity_and_source_outcomes.csv',C66/'all551_exact_reference_ORF_benchmark_models.csv',C66/'all553_independent_source_DNA_to_AA_identity_audit.csv',C34/'table13_new_eight_cell_fitness.csv',C34/'table14_structure_predictions.csv',R/'meta_clusters_2026-09-13/transcode_screen/phase1_Ribo-seq_ORFs.all.bed']
pins=[]
for p in sourcefiles+[REF]:
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 pins.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
assert pins[-1]['sha256']=='1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0'
(O/'PRE_COMPARISON_SOURCE_AND_CURRENT_REFERENCE_PINS.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pins':pins},indent=2))
old=pd.read_csv(sourcefiles[0],usecols=['source_Excel_row','Name','ORF type','ORF sequence','Amino Acid','Initiation codon','Ensembl/GENCODE Gene ID (v19)'],keep_default_na=False);models=pd.read_csv(sourcefiles[1],keep_default_na=False);audit=pd.read_csv(sourcefiles[2],keep_default_na=False);assert len(old)==553 and len(models)==551 and old.source_Excel_row.is_unique and models.source_Excel_row.is_unique
new=pd.read_csv(sourcefiles[3],usecols=['source_excel_row','cell_line','orf_targeted','type','biotype','Ensembl_id'],keep_default_na=False);new=new[new.type.eq('ncORF')];assert len(new)==17568 and new.orf_targeted.nunique()==2196
for col in ['Ensembl_id','biotype']:assert new.groupby('orf_targeted')[col].nunique().eq(1).all()
newtargets=new.sort_values('source_excel_row').drop_duplicates('orf_targeted');structure=pd.read_csv(sourcefiles[4],usecols=['source_excel_row','orf_id','sequence'],keep_default_na=False);assert structure.orf_id.is_unique
f=REF.open('rb');head=f.read(16);endian='<' if struct.unpack('<I',head[:4])[0]==0x1A412743 else '>';uint=lambda:struct.unpack(endian+'I',f.read(4))[0];sig,ver,count,res=struct.unpack(endian+'4I',head);assert sig==0x1A412743 and ver==0;idx={}
for _ in range(count):n=f.read(1)[0];chrom=f.read(n).decode();idx[chrom]=uint()
cache={}
def segment(chrom,a,b):
 if chrom not in cache:
  f.seek(idx[chrom]);size=uint();nn=uint();ns=[uint() for _ in range(nn)];nl=[uint() for _ in range(nn)];nm=uint();f.seek(8*nm,1);assert uint()==0;cache[chrom]=(size,list(zip(ns,nl)),f.tell())
 size,nb,offset=cache[chrom];assert 0<=a<b<=size;lo=a//4;f.seek(offset+lo);buf=np.frombuffer(f.read((b+3)//4-lo),dtype=np.uint8);p=np.arange(a,b);s=np.frombuffer(b'TCAG',dtype=np.uint8)[(buf[p//4-lo]>>(6-2*(p%4)))&3].copy()
 for start,length in nb:
  if start<b and start+length>a:s[max(0,start-a):min(b,start+length)-a]=ord('N')
 return s.tobytes().decode()
def reconstruct(chrom,blocks,strand):
 s=''.join(segment(chrom,a,b) for a,b in sorted(blocks));return s if strand=='+' else s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
code=dict(zip([''.join(t) for t in itertools.product('TCAG',repeat=3)],'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def translate(s):return ''.join(code.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def stripstop(s):return s[:-1] if s.endswith('*') else s
near={s for s in code if s!='ATG' and sum(a!=b for a,b in zip(s,'ATG'))==1}
def product_fields(dna,chrom,blocks,strand):
 chain=translate(dna);aa=stripstop(chain);inrna=sorted(blocks,reverse=strand=='-');start=inrna[0][0] if strand=='+' else inrna[0][1]-1
 return {'chrom':chrom,'strand':strand,'blocks_in_RNA_order0':'|'.join(str(a)+'-'+str(b) for a,b in inrna),'coding_first_base0':start,'DNA_sequence_in_source_stop_convention':dna,'DNA_nt':len(dna),'start_triplet':dna[:3],'terminal_inframe_stop':chain.endswith('*'),'internal_stop_count':aa.count('*'),'literal_standard_AA':aa,'AA_length':len(aa),'near_cognate_reference_start':dna[:3] in near,'initiator_M_sensitivity_AA':('M'+aa[1:]) if dna[:3] in near else aa,'RNA_chain_model_key':hashlib.sha256((chrom+':'+strand+':'+str(inrna)).encode()).hexdigest()}
def gene(s):return re.sub(r'\.\d+$','',s) if re.fullmatch(r'ENSG\d+(?:\.\d+)?',s) else ''
older=[];modelmap=models.set_index('source_Excel_row').to_dict('index');amap=audit.set_index('source_Excel_row').to_dict('index')
for r in old.to_dict('records'):
 rn=r['source_Excel_row'];sourceDNA=re.sub(r'\s+','',r['ORF sequence']).upper();sourceAA=re.sub(r'\s+','',r['Amino Acid']).upper();g=str(r['Ensembl/GENCODE Gene ID (v19)']);out={'source_Excel_row':rn,'source_name':r['Name'],'source_ORF_type':r['ORF type'],'source_gene_ID_literal':g,'source_single_stable_gene_ID':gene(g),'source_DNA_literal':sourceDNA,'source_printed_AA_literal':sourceAA,'source_printed_AA_without_one_terminal_stop':stripstop(sourceAA),'source_reported_start':r['Initiation codon'],'reference_mapped':rn in modelmap,'source_AA_other_difference_prior_flag':amap[rn]['source_AA_other_difference']}
 if rn in modelmap:
  m=modelmap[rn];blocks=[(int(a),int(b)) for a,b in re.findall(r'(\d+)-(\d+)',m['blocks_in_RNA_order0'])];dna=reconstruct(m['chrom'],blocks,m['strand']);assert dna==m['DNA_sequence']==sourceDNA;out.update(product_fields(dna,m['chrom'],blocks,m['strand']));out['older_benchmark_preserved_AA_literal']=m['AA_sequence'];out['literal_translation_equals_source_printed_AA']=stripstop(translate(dna))==stripstop(sourceAA);assert stripstop(m['AA_sequence'])==stripstop(sourceAA)
 older.append(out)
bed=collections.defaultdict(list)
for lineno,line in enumerate(sourcefiles[-1].read_text().splitlines(),1):
 v=line.split('\t')
 if len(v)<12:continue
 id=v[3].split(':')[0];sizes=list(map(int,v[10].rstrip(',').split(',')));starts=[int(v[1])+int(x) for x in v[11].rstrip(',').split(',')];blocks=list(zip(starts,[a+b for a,b in zip(starts,sizes)]));bed[id].append({'source_BED_line':lineno,'source_BED_name':v[3],'chrom':v[0],'strand':v[5],'blocks':blocks})
sm=structure.set_index('orf_id').to_dict('index');native=[]
for r in newtargets.to_dict('records'):
 id=r['orf_targeted'];out={'source_Excel_row':r['source_excel_row'],'source_ORF_ID':id,'source_biotype':r['biotype'],'source_gene_ID_literal':r['Ensembl_id'],'source_single_stable_gene_ID':gene(r['Ensembl_id']),'phase1_reference_model_count':len(bed.get(id,[])),'structure_reference_record_present':id in sm,'reference_mapped':len(bed.get(id,[]))==1}
 if out['reference_mapped']:
  m=bed[id][0];dna=reconstruct(m['chrom'],m['blocks'],m['strand']);out.update(product_fields(dna,m['chrom'],m['blocks'],m['strand']));out.update({k:v for k,v in m.items() if k not in ['blocks','chrom','strand']})
 if id in sm:out.update({'structure_source_Excel_row':sm[id]['source_excel_row'],'structure_full_AA_literal':sm[id]['sequence']})
 native.append(out)
pd.DataFrame(older).to_csv(O/'all553_original2021_complete_sequence_models_before_cross_screen_comparison.csv',index=False);pd.DataFrame(native).to_csv(O/'all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv',index=False)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'all_source_reference_models_prepared_before_new_identity_comparison','older_source_candidates':len(older),'older_reference_mapped':sum(x['reference_mapped'] for x in older),'older_unmapped_source_names':[x['source_name'] for x in older if not x['reference_mapped']],'older_single_literal_Ensembl_gene_IDs':sum(bool(x['source_single_stable_gene_ID']) for x in older),'newer_source_targets':len(native),'newer_unique_phase1_reference_models':sum(x['reference_mapped'] for x in native),'newer_unresolved_target_models':[{'id':x['source_ORF_ID'],'matches':x['phase1_reference_model_count']} for x in native if not x['reference_mapped']],'newer_single_literal_Ensembl_gene_IDs':sum(bool(x['source_single_stable_gene_ID']) for x in native),'newer_structure_reference_records':sum(x['structure_reference_record_present'] for x in native),'reference_N_characters_in_prepared_models':sum(x.get('DNA_sequence_in_source_stop_convention','').count('N') for x in older+native),'no_cross_screen_identity_or_outcome_comparisons_yet':True};(O/'COMPLETE_INPUT_REFERENCE_MODEL_PREFLIGHT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
