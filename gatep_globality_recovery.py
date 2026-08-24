#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, subprocess, time
from pathlib import Path
import yaml

LCDM={"logA":3.0516831,"n_s":0.97006392,"H0":67.003038,"omega_b":0.022462565,"omega_cdm":0.12107093,"tau_reio":0.056456949,"Tcal":0.99859358,"Ecal":1.0066561,"A_act":1.0011501,"P_act":1.0016445}
OLD1={"logA":3.0503134,"n_s":0.97792072,"H0":68.510891,"omega_b":0.022477741,"omega_cdm":0.12930926,"tau_reio":0.050513103,"Tcal":0.99898779,"Ecal":1.0047459,"log10_maxion_ac":-3.5342913,"fraction_maxion_ac":0.065121699,"A_act":1.0003491,"P_act":1.0004015}
OLD2={"logA":3.0496313,"n_s":0.97943114,"H0":68.982829,"omega_b":0.022480985,"omega_cdm":0.12977741,"tau_reio":0.051216615,"Tcal":0.99998665,"Ecal":1.0027784,"log10_ac1":-3.4759976,"log10_ac2":-4.0221175,"f_ax1":0.056720487,"f_ax2":0.025442727,"A_act":1.000423,"P_act":0.99971618}
AC1=(-3.7500772228,-3.0004340775); AC2=(-4.5000137334,-3.7500772228); ACA=(-4.5000137334,-3.0004340775)
def mid(x): return (x[0]+x[1])/2

def base(model):
    extra={"non linear":"hmcode","N_ur":2.0328,"N_ncdm":1,"T_ncdm":0.71611,"lensing":"yes","P_k_max_h/Mpc":1.0,"do_shooting":"yes","do_shooting_mscf":"yes","attractor_ic_scf":"no","loop_over_background_for_closure_relation":"no","background_Nloga":100000,"tol_shooting_deltaF":0.01,"tol_shooting_deltax":0.01}
    if model=='lcdm': extra['N_mscf']=0
    elif model=='1ax': extra.update({"N_mscf":1,"n_axion_mscf":"3","theta_ini_mscf":"2.8","theta_prime_ini_mscf":"0."})
    elif model=='2ax': extra.update({"N_mscf":2,"n_axion_mscf":"3,3","theta_ini_mscf":"2.8,2.8","theta_prime_ini_mscf":"0.,0."})
    p={
      "logA":{"prior":{"min":1.61,"max":3.91},"ref":3.05,"proposal":0.002,"drop":True},
      "A_s":{"value":"lambda logA: 1e-10*np.exp(logA)"},
      "n_s":{"prior":{"min":0.8,"max":1.2},"ref":0.98,"proposal":0.003},
      "H0":{"prior":{"min":20.0,"max":100.0},"ref":67.0,"proposal":1.0},
      "omega_b":{"prior":{"min":0.005,"max":0.1},"ref":0.02246,"proposal":0.0001},
      "omega_cdm":{"prior":{"min":0.001,"max":0.99},"ref":0.121,"proposal":0.001},
      "m_ncdm":{"value":0.06},
      "tau_reio":{"prior":{"dist":"norm","loc":0.051,"scale":0.006},"ref":0.051,"proposal":0.003},
      "Tcal":{"prior":{"min":0.8,"max":1.2},"ref":1.0,"proposal":0.003},
      "Ecal":{"prior":{"min":0.8,"max":1.2},"ref":1.0,"proposal":0.003}}
    if model=='1ax':
      p.update({"log10_maxion_ac":{"prior":{"min":ACA[0],"max":ACA[1]},"ref":mid(ACA),"proposal":0.08},"fraction_maxion_ac":{"prior":{"min":0.0,"max":0.3},"ref":0.05,"proposal":0.015}})
    if model=='2ax':
      p.update({
       "log10_ac1":{"prior":{"min":AC1[0],"max":AC1[1]},"ref":mid(AC1),"proposal":0.06,"drop":True},
       "log10_ac2":{"prior":{"min":AC2[0],"max":AC2[1]},"ref":mid(AC2),"proposal":0.08,"drop":True},
       "f_ax1":{"prior":{"min":0.0,"max":0.3},"ref":0.05,"proposal":0.015,"drop":True},
       "f_ax2":{"prior":{"min":0.0,"max":0.3},"ref":0.025,"proposal":0.015,"drop":True},
       "log10_maxion_ac":{"value":'lambda log10_ac1, log10_ac2: "%.16g,%.16g" % (log10_ac1, log10_ac2)',"derived":False},
       "fraction_maxion_ac":{"value":'lambda f_ax1, f_ax2: "%.16g,%.16g" % (f_ax1 if f_ax1 > 1e-12 else 1e-12, f_ax2 if f_ax2 > 1e-12 else 1e-12)',"derived":False}})
    return {
      "theory":{"classy":{"path":str(Path.cwd()/"external/mAxiCLASS"),"ignore_obsolete":True,"extra_args":extra}},
      "likelihood":{
       "act_dr6_cmbonly.ACTDR6CMBonly":{"input_file":"dr6_data_cmbonly.fits","lmax_theory":9000,"ell_cuts":{"TT":[600,8500],"TE":[600,8500],"EE":[600,8500]},"stop_at_error":True,"params":{"A_act":{"prior":{"min":0.5,"max":1.5},"ref":1.0,"proposal":0.003},"P_act":{"prior":{"min":0.9,"max":1.1},"ref":1.0,"proposal":0.01}}},
       "candl_like":{"external":"__CANDL_EXTERNAL__","data_set_file":str(Path.cwd()/"external/spt_candl_data/spt_candl_data/SPT3G_D1_TnE_v0/SPT3G_D1_TnE_index.yaml"),"variant":"lite","clear_internal_priors":True,"additional_args":{},"feedback":True,"wrapper":None}},
      "prior":{"cal_dip_prior":"lambda A_act: stats.norm.logpdf(A_act, loc=1.0, scale=0.003)","gaussian_Tcal":"lambda Tcal: stats.norm.logpdf(Tcal, loc=1.0, scale=0.0036)"},
      "params":p,"packages_path":str(Path.cwd()/".cobaya"),"timing":True}

def set_refs(info,r):
    for k,v in r.items():
      if k in ('A_act','P_act'): info['likelihood']['act_dr6_cmbonly.ACTDR6CMBonly']['params'][k]['ref']=float(v)
      elif k in info['params'] and isinstance(info['params'][k],dict) and 'prior' in info['params'][k]: info['params'][k]['ref']=float(v)

def dump_yaml(path,info):
    t=yaml.safe_dump(info,sort_keys=False,width=120)
    t=t.replace("external: __CANDL_EXTERNAL__","external: !!python/name:candl.interface.CandlCobayaLikelihood ''")
    t=t.replace("external: '__CANDL_EXTERNAL__'","external: !!python/name:candl.interface.CandlCobayaLikelihood ''")
    path.write_text(t)

def parse_table(path):
    lines=[x.strip() for x in path.read_text().splitlines() if x.strip()]
    h=lines[0].lstrip('#').split(); v=lines[1].split()
    if len(h)!=len(v): raise RuntimeError(f'column mismatch {path}')
    out={}
    for a,b in zip(h,v):
      try: out[a]=float(b)
      except ValueError: pass
    return out

def effective(row):
    pen={"tau":((row['tau_reio']-0.051)/0.006)**2,"A_act_cal":((row['A_act']-1)/0.003)**2,"Tcal":((row['Tcal']-1)/0.0036)**2}
    return {"likelihood_chi2":row['chi2'],"gaussian_penalties":pen,"effective_chi2":row['chi2']+sum(pen.values()),"chi2_components":{k:v for k,v in row.items() if k.startswith('chi2__')},"parameters":{k:v for k,v in row.items() if not k.startswith('minuslog') and not k.startswith('chi2') and k!='weight'}}

def run_candidate(model,label,r,wd,kind,max_evals=650,rhoend=0.03):
    d=wd/'runs'/model/label; d.mkdir(parents=True,exist_ok=True)
    info=base(model); set_refs(info,r); info['output']=str(d/label)
    if kind=='evaluate': info['sampler']={'evaluate':{}}
    else: info['sampler']={'minimize':{'method':'bobyqa','ignore_prior':False,'max_evals':max_evals,'best_of':1,'seed':sum(map(ord,model+label))*7919 % 2000000000,'override_bobyqa':{'rhoend':rhoend}}}
    y=d/f'{label}.yaml'; dump_yaml(y,info); log=d/f'{label}.log'; t=time.time()
    with log.open('w') as f: p=subprocess.run(['cobaya-run',str(y)],stdout=f,stderr=subprocess.STDOUT,text=True)
    out={"model":model,"label":label,"kind":kind,"refs":r,"exit_code":p.returncode,"runtime_seconds":time.time()-t,"yaml":str(y),"log":str(log),"valid":False}
    if p.returncode: return out
    try:
      txt=d/(f'{label}.1.txt' if kind=='evaluate' else f'{label}.minimum.txt')
      out.update(effective(parse_table(txt))); out['result_file']=str(txt); out['valid']=True
    except Exception as e: out['parse_error']=repr(e)
    return out

def save(wd,m,c,b):
    o=wd/'results'/m; o.mkdir(parents=True,exist_ok=True)
    (o/'candidates.json').write_text(json.dumps(c,indent=2)); (o/'best.json').write_text(json.dumps(b,indent=2))
    with (o/'candidates.csv').open('w',newline='') as f:
      keys=['label','kind','valid','exit_code','runtime_seconds','effective_chi2','likelihood_chi2']; w=csv.DictWriter(f,fieldnames=keys); w.writeheader()
      for x in c: w.writerow({k:x.get(k) for k in keys})

def refine(model,c,wd):
    mins=[x for x in c if x.get('valid') and x['kind']=='minimize']
    if mins:
      s=min(mins,key=lambda x:x['effective_chi2'])
      c.append(run_candidate(model,'refine_best_minimum',s['parameters'],wd,'minimize',1100,0.005))
    return min([x for x in c if x.get('valid')],key=lambda x:x['effective_chi2'])

def stage1(wd):
    lc=[run_candidate('lcdm','lcdm_anchor_eval',LCDM,wd,'evaluate')]; lv=[x for x in lc if x.get('valid')]
    if not lv: raise SystemExit('LCDM anchor failed')
    lb=min(lv,key=lambda x:x['effective_chi2']); save(wd,'lcdm',lc,lb)
    c=[]
    for i,a in enumerate((ACA[0]+0.05,mid(ACA),ACA[1]-0.05),1): c.append(run_candidate('1ax',f'boundary_eval_{i}',{**LCDM,'log10_maxion_ac':a,'fraction_maxion_ac':0.0},wd,'evaluate'))
    c.append(run_candidate('1ax','near_boundary_min',{**LCDM,'log10_maxion_ac':mid(ACA),'fraction_maxion_ac':1e-8},wd,'minimize'))
    c.append(run_candidate('1ax','old_local_min',OLD1,wd,'minimize'))
    c.append(run_candidate('1ax','interior_early_min',{**LCDM,'H0':69.0,'omega_cdm':0.128,'log10_maxion_ac':-3.20,'fraction_maxion_ac':0.05},wd,'minimize'))
    c.append(run_candidate('1ax','interior_late_min',{**LCDM,'H0':69.0,'omega_cdm':0.128,'log10_maxion_ac':-4.25,'fraction_maxion_ac':0.05},wd,'minimize'))
    if not [x for x in c if x.get('valid')]: raise SystemExit('No valid 1ax candidate')
    b=refine('1ax',c,wd); save(wd,'1ax',c,b)
    ok=b['effective_chi2']<=lb['effective_chi2']+0.1
    (wd/'results/stage1_status.json').write_text(json.dumps({'lcdm':lb['effective_chi2'],'1ax':b['effective_chi2'],'nesting_1ax_le_lcdm_tol0p1':ok},indent=2))
    if not ok: raise SystemExit('1ax nesting recovery failed')

def embed(one):
    p=one['parameters']; r={k:p[k] for k in LCDM if k in p}; a=p.get('log10_maxion_ac',mid(AC1)); f=p.get('fraction_maxion_ac',0.0)
    if AC1[0]<=a<=AC1[1]: r.update({'log10_ac1':a,'f_ax1':f,'log10_ac2':mid(AC2),'f_ax2':0.0})
    elif AC2[0]<=a<=AC2[1]: r.update({'log10_ac1':mid(AC1),'f_ax1':0.0,'log10_ac2':a,'f_ax2':f})
    else: raise RuntimeError('1ax ac outside frozen 2ax union')
    return r

def stage2(wd):
    lb=json.loads((wd/'results/lcdm/best.json').read_text()); ob=json.loads((wd/'results/1ax/best.json').read_text()); c=[]
    for i,(a1,a2) in enumerate(((AC1[0]+0.05,AC2[0]+0.05),(mid(AC1),mid(AC2)),(AC1[1]-0.05,AC2[1]-0.05)),1): c.append(run_candidate('2ax',f'lcdm_boundary_eval_{i}',{**LCDM,'log10_ac1':a1,'log10_ac2':a2,'f_ax1':0.0,'f_ax2':0.0},wd,'evaluate'))
    e=embed(ob); c.append(run_candidate('2ax','oneax_embedding_eval',e,wd,'evaluate'))
    en=dict(e); en['f_ax1']=max(en['f_ax1'],1e-8); en['f_ax2']=max(en['f_ax2'],1e-8)
    c.append(run_candidate('2ax','oneax_embedding_min',en,wd,'minimize'))
    c.append(run_candidate('2ax','old_local_min',OLD2,wd,'minimize'))
    c.append(run_candidate('2ax','asym_f1_dominant_min',{**LCDM,'H0':69.0,'omega_cdm':0.129,'log10_ac1':-3.35,'log10_ac2':-4.05,'f_ax1':0.08,'f_ax2':0.01},wd,'minimize'))
    c.append(run_candidate('2ax','asym_f2_dominant_min',{**LCDM,'H0':69.0,'omega_cdm':0.129,'log10_ac1':-3.55,'log10_ac2':-4.15,'f_ax1':0.02,'f_ax2':0.08},wd,'minimize'))
    if not [x for x in c if x.get('valid')]: raise SystemExit('No valid 2ax candidate')
    b=refine('2ax',c,wd); save(wd,'2ax',c,b)
    ok1=ob['effective_chi2']<=lb['effective_chi2']+0.1; ok2=b['effective_chi2']<=ob['effective_chi2']+0.1
    (wd/'results/stage2_status.json').write_text(json.dumps({'lcdm':lb['effective_chi2'],'1ax':ob['effective_chi2'],'2ax':b['effective_chi2'],'nesting_1ax_le_lcdm_tol0p1':ok1,'nesting_2ax_le_1ax_tol0p1':ok2},indent=2))
    if not(ok1 and ok2): raise SystemExit('2ax nesting recovery failed')

def report(wd):
    l=json.loads((wd/'results/lcdm/best.json').read_text()); o=json.loads((wd/'results/1ax/best.json').read_text()); t=json.loads((wd/'results/2ax/best.json').read_text())
    cl,co,ct=l['effective_chi2'],o['effective_chi2'],t['effective_chi2']; d=ct-co; a=d+4; f2=t['parameters'].get('f_ax2'); fb=f2 is not None and f2<=1e-5; nesting=ct<=co+0.1 and co<=cl+0.1
    route='NUMERICAL WARNING — nesting failed; no Gate-P verdict.' if not nesting else ('PROVISIONAL POSITIVE GATE P — run theory-precision gate, then Gate F.' if d<-4 and not fb else 'PROVISIONAL NEGATIVE/NO-COMPLEXITY GATE P — run theory-precision gate, then Virkelighedstest.')
    r={'case':'CASE-Q013','stage':'GATE-P GLOBALITY RECOVERY','effective_chi2':{'lcdm':cl,'1ax':co,'2ax':ct},'nesting':{'2ax_le_1ax':ct<=co+0.1,'1ax_le_lcdm':co<=cl+0.1,'tolerance':0.1},'delta_chi2_P_2ax_minus_1ax':d,'delta_AIC_P':a,'f2':f2,'f2_boundary_le_1e-5':fb,'theory_precision_gate':'NOT_RUN_BY_THIS_WORKFLOW','route':route,'important':'No Wilks/sigma conversion. Final Gate-P verdict requires |delta chi2_theory| <= 0.1.'}
    (wd/'GATEP_GLOBALITY_RESULT.json').write_text(json.dumps(r,indent=2))
    md=(f"# CASE-Q013 — Gate-P Globality Recovery\n\n"
        f"- LCDM effective chi2: **{cl:.9f}**\n"
        f"- 1ax effective chi2: **{co:.9f}**\n"
        f"- 2ax effective chi2: **{ct:.9f}**\n"
        f"- Nesting 1ax <= LCDM (+0.1): **{r['nesting']['1ax_le_lcdm']}**\n"
        f"- Nesting 2ax <= 1ax (+0.1): **{r['nesting']['2ax_le_1ax']}**\n"
        f"- Delta chi2 P: **{d:+.9f}**\n"
        f"- Delta AIC P: **{a:+.9f}**\n"
        f"- f2: **{f2}**\n"
        f"- f2 boundary (<=1e-5): **{fb}**\n"
        f"- Theory precision gate: **NOT RUN**\n"
        f"- Route: **{route}**\n")
    (wd/'GATEP_GLOBALITY_REPORT.md').write_text(md)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--stage',choices=['1ax','2ax','report'],required=True); ap.add_argument('--workdir',default='gatep_recovery'); a=ap.parse_args(); wd=Path(a.workdir); (wd/'results').mkdir(parents=True,exist_ok=True)
    {'1ax':stage1,'2ax':stage2,'report':report}[a.stage](wd)
if __name__=='__main__': main()
