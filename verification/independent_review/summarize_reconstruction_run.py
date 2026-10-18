"""Collect the active replay jobs, preserving failed experimental run records."""
from pathlib import Path
from datetime import datetime,timezone
import re,json,hashlib


def main():
    base=Path(__file__).parent;root=base.parent.parent
    wrapper=(base/'run_independent_review.ps1').read_text(encoding='utf-8')
    jobs=re.findall(r"^\s*'([^']+\.py)'",wrapper,re.MULTILINE)
    records=[]
    for name in jobs:
        p=base/'run_logs'/Path(name).with_suffix('.run.json')
        data=json.loads(p.read_text(encoding='utf-8'))
        records.append({'script':name,'source_sha256':hashlib.sha256((base/name).read_bytes()).hexdigest(),
                        'exit_code':data['exit_code'],'elapsed_seconds':data['elapsed_seconds'],
                        'peak_rss_mb':data['peak_rss_mb'],
                        'memory_limit_mb':data['memory_limit_mb'],
                        'windows_job_memory_limit_enabled':data['windows_job_memory_limit_enabled'],
                        'terminated_reason':data['terminated_reason']})
    files=[p for p in base.glob('*.json') if p.name not in ('reconstruction_run_summary.json',)]
    files += [root/'main_merged_revised.tex',root/'main_merged_revised.pdf',
              root/'remaining_proof_obligations_zh.tex',root/'remaining_proof_obligations_zh.pdf']
    files += [base.parent/'requirements.lock.txt',base.parent/'environment.json',base/'run_certificate_job.py',base/'run_independent_review.ps1',root/'compile_revised_pdfs.ps1']
    hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),
            'active_jobs':len(jobs),'successful_active_jobs':sum(r['exit_code']==0 for r in records),
            'all_active_jobs_successful':all(r['exit_code']==0 for r in records),
            'sum_worker_seconds':round(sum(r['elapsed_seconds'] for r in records),3),
            'maximum_worker_peak_rss_mb':max(r['peak_rss_mb'] for r in records),
            'jobs':records,'artifact_sha256':hashes,
            'scope':'Execution success applies to the listed implemented reconstructions, not to every mathematical assertion. See audit report for remaining proof obligations.'}
    out=base/'reconstruction_run_summary.json'
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('active_jobs','successful_active_jobs','all_active_jobs_successful','sum_worker_seconds','maximum_worker_peak_rss_mb')},indent=2))


if __name__=='__main__':main()
