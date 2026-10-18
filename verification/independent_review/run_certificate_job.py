"""Run one certificate worker with a memory limit and persistent JSON log."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
import psutil


def windows_memory_job(worker,limit_bytes):
    """Put only this new worker and its children in a Windows memory job."""
    if os.name!='nt':return None
    import ctypes as c
    from ctypes import wintypes as w
    class Basic(c.Structure):
        _fields_=[('PerProcessUserTimeLimit',c.c_longlong),('PerJobUserTimeLimit',c.c_longlong),
                  ('LimitFlags',w.DWORD),('MinimumWorkingSetSize',c.c_size_t),
                  ('MaximumWorkingSetSize',c.c_size_t),('ActiveProcessLimit',w.DWORD),
                  ('Affinity',c.c_size_t),('PriorityClass',w.DWORD),('SchedulingClass',w.DWORD)]
    class IO(c.Structure):
        _fields_=[(k,c.c_ulonglong) for k in ('ReadOperationCount','WriteOperationCount',
                  'OtherOperationCount','ReadTransferCount','WriteTransferCount','OtherTransferCount')]
    class Extended(c.Structure):
        _fields_=[('BasicLimitInformation',Basic),('IoInfo',IO),
                  ('ProcessMemoryLimit',c.c_size_t),('JobMemoryLimit',c.c_size_t),
                  ('PeakProcessMemoryUsed',c.c_size_t),('PeakJobMemoryUsed',c.c_size_t)]
    kernel=c.WinDLL('kernel32',use_last_error=True)
    kernel.CreateJobObjectW.argtypes=[c.c_void_p,w.LPCWSTR];kernel.CreateJobObjectW.restype=w.HANDLE
    kernel.SetInformationJobObject.argtypes=[w.HANDLE,c.c_int,c.c_void_p,w.DWORD]
    kernel.AssignProcessToJobObject.argtypes=[w.HANDLE,w.HANDLE]
    kernel.CloseHandle.argtypes=[w.HANDLE]
    job=kernel.CreateJobObjectW(None,None)
    if not job:return None
    info=Extended();info.BasicLimitInformation.LimitFlags=0x200|0x2000
    info.JobMemoryLimit=limit_bytes
    if not kernel.SetInformationJobObject(job,9,c.byref(info),c.sizeof(info)) or not kernel.AssignProcessToJobObject(job,int(worker._handle)):
        kernel.CloseHandle(job);return None
    return job


def main():
    p=argparse.ArgumentParser()
    p.add_argument('script');p.add_argument('--memory-mb',type=int,default=1024)
    p.add_argument('--timeout-s',type=int,default=900)
    a=p.parse_args()
    script=Path(a.script).resolve();logs=script.parent/'run_logs';logs.mkdir(exist_ok=True)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    start=time.monotonic();peak=0;reason=None
    output=logs/(script.stem+'.log')
    with output.open('w',encoding='utf-8') as f:
        worker=subprocess.Popen([sys.executable,'-u',str(script)],cwd=script.parent,env=env,
                                stdout=f,stderr=subprocess.STDOUT)
        memory_job=windows_memory_job(worker,a.memory_mb*1024**2)
        process=psutil.Process(worker.pid)
        while worker.poll() is None:
            try:
                rss=process.memory_info().rss
                for child in process.children(recursive=True):rss+=child.memory_info().rss
                peak=max(peak,rss)
                if rss>a.memory_mb*1024**2:reason='memory limit exceeded'
                if time.monotonic()-start>a.timeout_s:reason='time limit exceeded'
                if reason:
                    for child in process.children(recursive=True):child.terminate()
                    worker.terminate();break
            except psutil.NoSuchProcess:break
            time.sleep(0.25)
        rc=worker.wait()
        if memory_job:
            import ctypes
            ctypes.windll.kernel32.CloseHandle(ctypes.c_void_p(memory_job))
    result={'script':script.name,'python':sys.executable,'exit_code':rc,
            'elapsed_seconds':round(time.monotonic()-start,3),
            'peak_rss_mb':round(peak/1024**2,3),'memory_limit_mb':a.memory_mb,
            'windows_job_memory_limit_enabled':bool(memory_job),
            'timeout_s':a.timeout_s,'terminated_reason':reason,'stdout_log':str(output)}
    (logs/(script.stem+'.run.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    print(output.read_text(encoding='utf-8')[-4000:])
    sys.exit(rc if rc else (1 if reason else 0))


if __name__=='__main__':main()
