from pathlib import Path
import argparse,datetime,hashlib,json,re
import numpy as np
import pandas as pd
M=Path(__file__).resolve().parent
N=M.parent if (M.parent/'cycle_87_fresh_guide_cross_screen_nORF_validation').exists() else M.parents[1]/'exploration_2026-10-11'/'novelty_discovery'
C87=N/'cycle_87_fresh_guide_cross_screen_nORF_validation'
CELLS=['A375','HEPG2']
SOURCE=M/'sources'/'41587_2020_806_MOESM5_ESM.xlsx'
COUNTS=M/'sources'/'41587_2020_806_MOESM12_ESM.xlsx'
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def savej(name,value):
 (M/name).write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf8')
def pin(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def flag(s):return s.astype(str).str.lower().eq('true')
def prepare():
 meta=pd.read_excel(SOURCE,sheet_name='Table S27',header=2,usecols=['Construct Barcode','Gene Symbol','Gene_type','Type'])
 lib=pd.read_excel(SOURCE,sheet_name='Table S24',header=2)
 raw=pd.read_excel(COUNTS)
 assert raw['Construct Barcode'].is_unique
 assert meta['Construct Barcode'].is_unique
 coverage_cols=[x for x in raw if 'ETP' in x and (x.startswith('A375 ') or x.startswith('HepG2 '))]
 assert len(coverage_cols)==6
 rawmeta=raw[['Construct Barcode','pDNA']+coverage_cols].copy()
 rawmeta['secondary_ETP_all6_ge10']=rawmeta[coverage_cols].ge(10).all(axis=1)
 rawmeta=rawmeta[['Construct Barcode','pDNA','secondary_ETP_all6_ge10']]
 lib_counts=lib.groupby('sgRNA Sequence').agg(secondary_library_rows=('Input','size'),secondary_library_assignment_count=('Input','nunique')).reset_index().rename(columns={'sgRNA Sequence':'sgrna'})
 newer=pd.read_csv(C87/'all720_newer_exact_product_guide_roles_and_eligibility_before_outcome_join.csv')
 groups=pd.read_csv(C87/'all127_complete_product_outcome_eligibility_and_denominators.csv')[['pair_ID','dependence_cluster']]
 assert newer[['pair_ID','sgrna']].duplicated().sum()==0
 a=newer.merge(groups,on='pair_ID',validate='many_to_one').merge(meta,left_on='sgrna',right_on='Construct Barcode',how='left',validate='many_to_one').merge(rawmeta,left_on='sgrna',right_on='Construct Barcode',how='left',validate='many_to_one',suffixes=('','_counts')).merge(lib_counts,on='sgrna',how='left',validate='many_to_one')
 a['secondary_base_name']=a['Gene Symbol'].fillna('').str.replace(r'_exon[0-9]+$','',regex=True)
 a['literal_source_name_match']=a['secondary_base_name'].eq(a['old_source_name'])
 a['secondary_metadata_eligible']=flag(a['primary_eligible']) & a['Gene_type'].eq('ORF_exon') & a['literal_source_name_match'] & a['pDNA'].ge(10) & a['secondary_library_assignment_count'].eq(1)
 a.to_csv(M/'all720_followup_guide_metadata_eligibility.csv',index=False)
 universe=groups.merge(newer[['pair_ID','old_source_name','new_source_ORF_ID']].drop_duplicates(),on='pair_ID',validate='one_to_one')
 cnt=a[a.secondary_metadata_eligible].groupby('pair_ID').size()
 universe['eligible_shared_guides']=universe.pair_ID.map(cnt).fillna(0).astype(int)
 universe['min3_metadata_eligible']=universe.eligible_shared_guides.ge(3)
 universe.to_csv(M/'all127_followup_product_metadata_cohorts.csv',index=False)
 pins=[SOURCE,COUNTS,C87/'all720_newer_exact_product_guide_roles_and_eligibility_before_outcome_join.csv',C87/'all127_complete_product_outcome_eligibility_and_denominators.csv',C87/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv',M/'FROZEN_FOLLOWUP_GUIDE_PROTOCOL.json']
 savej('METADATA_PREFLIGHT.json',{'utc':utc(),'outcomes_joined':False,'source_count_rows':len(raw),'source_library_rows':len(lib),'unique_source_counts':int(raw['Construct Barcode'].nunique()),'secondary_exact_sequence_matches':int(a['Gene Symbol'].notna().sum()),'eligible_shared_guides_before_min3':int(a.secondary_metadata_eligible.sum()),'eligible_products_min3':int(universe.min3_metadata_eligible.sum()),'eligible_shared_guides_min3':int(universe.loc[universe.min3_metadata_eligible,'eligible_shared_guides'].sum()),'metadata_pins':[pin(p) for p in pins],'source_name_mismatch_rows':a[a['Gene Symbol'].notna() & ~a.literal_source_name_match][['sgrna','old_source_name','Gene Symbol']].to_dict('records')})
 print((M/'METADATA_PREFLIGHT.json').read_text())
def residual(d,xcol='secondary_LFC',ycol='newer_LFC',min_guides=3):
 d=d.copy()
 assert not d[['pair_ID','sgrna','cell']].duplicated().any()
 n=d.groupby('pair_ID').sgrna.nunique()
 ids=n[n>=min_guides].index
 d=d[d.pair_ID.isin(ids)].copy()
 counts=d.groupby(['pair_ID','sgrna']).size()
 assert counts.eq(d.cell.nunique()).all()
 for src,col in [('x',xcol),('y',ycol)]:d[src]=d[col]-d.groupby(['pair_ID','cell'])[col].transform('mean')
 d['weight']=1/(d.cell.nunique()*d.pair_ID.map(n))
 return d

def corr(d):
 w=d.weight.to_numpy();x=d.x.to_numpy();y=d.y.to_numpy()
 return float(np.sum(w*x*y)/np.sqrt(np.sum(w*x*x)*np.sum(w*y*y)))
def moments(d):
 z=d[['pair_ID','dependence_cluster']].copy()
 z['xx']=d.weight*d.x*d.x;z['yy']=d.weight*d.y*d.y;z['xy']=d.weight*d.x*d.y
 return z.groupby('dependence_cluster')[['xx','yy','xy']].sum()
def describe(d,label):return {'analysis':label,'products':int(d.pair_ID.nunique()),'guides':int(d.sgrna.nunique()),'guide_cell_pairs':len(d),'correlation':corr(d)}
def analyze():
 assert (M/'METADATA_PREFLIGHT.json').exists()
 a=pd.read_csv(M/'all720_followup_guide_metadata_eligibility.csv')
 eligible=a[flag(a.secondary_metadata_eligible)].copy()
 s=pd.read_excel(SOURCE,sheet_name='Table S27',header=2)
 raw=pd.read_excel(COUNTS)
 assert np.isfinite(raw.select_dtypes('number')).all().all()
 count_cols=list(raw.columns[2:])
 cpm=raw[count_cols]/raw[count_cols].sum(axis=0)*1e6
 logs=np.log2(cpm+1)
 source_rows=[]
 for cell,prefix,pub,etp in [('A375','A375 ','LFC A375 pDNA_D21','LFC A375 21D ETP'),('HEPG2','HepG2 ','LFC HepG2 pDNA_D21','LFC HepG2 21D ETP')]:
  end=[c for c in raw if c.startswith(prefix) and c.endswith('D21')]
  assert len(end)==3
  calc=logs[end].mean(axis=1)-logs['pDNA']
  r=s[['Construct Barcode','Construct IDs',pub,etp]].rename(columns={'Construct Barcode':'sgrna','Construct IDs':'secondary_construct_id',pub:'secondary_LFC',etp:'secondary_ETP_LFC'})
  r['cell']=cell
  direct=pd.DataFrame({'sgrna':raw['Construct Barcode'],'secondary_direct_meanlog_CPM_LFC':calc})
  r=r.merge(direct,on='sgrna',validate='one_to_one')
  source_rows.append(r)
 src=pd.concat(source_rows,ignore_index=True)
 src.to_csv(M/'all13992_secondary_source_LFCs_and_direct_count_reconstruction.csv',index=False)
 newer=pd.read_csv(C87/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv')
 newer=newer[newer.cell.isin(CELLS)][['pair_ID','sgrna','cell','LFC']].rename(columns={'LFC':'newer_LFC'})
 d=eligible.merge(src,on='sgrna',validate='one_to_many').merge(newer,on=['pair_ID','sgrna','cell'],validate='one_to_one')
 finite=np.isfinite(d[['secondary_LFC','newer_LFC']]).all(axis=1)
 guidefinite=d.assign(finite=finite).groupby(['pair_ID','sgrna']).finite.all()
 invalid=set(guidefinite[~guidefinite].index)
 d['finite_in_both_cells']=[(p,g) not in invalid for p,g in zip(d.pair_ID,d.sgrna)]
 d.to_csv(M/'all_eligible_followup_outcomes_before_min3.csv',index=False)
 d=d[d.finite_in_both_cells].copy()
 main=residual(d)
 assert len(main)>0
 main.to_csv(M/'primary_followup_paired_residuals.csv',index=False)
 primary=describe(main,'primary_source_pDNA_to_D21')
 mo=moments(main);v=mo.to_numpy();k=len(v)
 rng=np.random.default_rng(20261016);draws=rng.integers(0,k,size=(2000,k))
 boot=v[draws].sum(axis=1);br=boot[:,2]/np.sqrt(boot[:,0]*boot[:,1])
 np.savez_compressed(M/'all2000_followup_cluster_draws.npz',draws=draws,cluster_ids=mo.index.to_numpy(dtype=str))
 pd.DataFrame({'draw':np.arange(2000),'correlation':br}).to_csv(M/'all2000_followup_cluster_bootstrap.csv',index=False)
 primary['dependence_clusters']=k;primary['nominal95_cluster_CI']=np.quantile(br,[.025,.975]).tolist()
 arrays=[]
 for pid,q in main.groupby('pair_ID',sort=True):
  x=q.pivot(index='sgrna',columns='cell',values='x').reindex(columns=CELLS).to_numpy()
  y=q.pivot(index='sgrna',columns='cell',values='y').reindex(columns=CELLS).to_numpy()
  arrays.append((x,y,1/(len(x)*len(CELLS))))
 denom=np.sqrt(float((main.weight*main.x**2).sum()*(main.weight*main.y**2).sum()))
 rng=np.random.default_rng(20261017);perms=np.zeros(5000)
 for j in range(5000):perms[j]=sum(float((x*y[rng.permutation(len(y))]).sum())*w for x,y,w in arrays)/denom
 pd.DataFrame({'draw':np.arange(5000),'correlation':perms}).to_csv(M/'all5000_followup_identity_permutations.csv',index=False)
 primary['permutation_extreme_count']=int(np.sum(np.abs(perms)>=abs(primary['correlation'])))
 primary['permutation_two_sided_add_one_P']=(1+primary['permutation_extreme_count'])/5001
 sensitivities=[primary.copy()]
 for label,q,xcol,nmin in [('atleast4_shared',d,'secondary_LFC',4),('all6_secondary_ETP_ge10',d[flag(d.secondary_ETP_all6_ge10)],'secondary_LFC',3),('ETP_to_D21_source',main,'secondary_ETP_LFC',3),('direct_count_meanlogCPM',main,'secondary_direct_meanlog_CPM_LFC',3)]:
  sub=residual(q,xcol=xcol,min_guides=nmin)
  sensitivities.append(describe(sub,label))
 for cell in CELLS:sensitivities.append(describe(residual(main[main.cell.eq(cell)]),'cell_'+cell))
 omissions=[]
 for pid in main.pair_ID.unique():omissions.append({'omitted_pair_ID':pid,'correlation':corr(main[main.pair_ID.ne(pid)])})
 pd.DataFrame(omissions).to_csv(M/'all_single_product_omissions.csv',index=False)
 pd.DataFrame(sensitivities).to_csv(M/'all_declared_followup_sensitivities.csv',index=False)
 summary={'utc':utc(),'status':'analyzed_pending_separate_verification','primary':primary,'sensitivities':sensitivities,'single_product_omission_range':[min(q['correlation'] for q in omissions),max(q['correlation'] for q in omissions)],'source_direct_count_LFC_comparison':{'max_absolute_difference':float(np.max(np.abs(main.secondary_LFC-main.secondary_direct_meanlog_CPM_LFC))),'Pearson':float(main[['secondary_LFC','secondary_direct_meanlog_CPM_LFC']].corr().iloc[0,1])},'limitations':json.loads((M/'FROZEN_FOLLOWUP_GUIDE_PROTOCOL.json').read_text())['limits'],'software':{'numpy':np.__version__,'pandas':pd.__version__}}
 savej('FOLLOWUP_GUIDE_RESULTS.json',summary)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','analyze']);args=p.parse_args()
 (prepare if args.stage=='prepare' else analyze)()
