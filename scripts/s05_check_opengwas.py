import json, os, sys
rep = []
p = r'D:\endometriosis_project\.ieugwaspy.json'
rep.append(f"jwt file exists: {os.path.exists(p)}")
if os.path.exists(p):
    j = json.load(open(p))
    rep.append(f"keys: {list(j.keys())}")
    # 不打印 token 全文
    for k, v in j.items():
        if isinstance(v, str):
            rep.append(f"  {k}: len={len(v)} head={v[:12]}...")
        else:
            rep.append(f"  {k}: {v}")
try:
    import ieugwasr as ig
    rep.append(f"ieugwaspy version: {getattr(ig,'__version__','?')}")
except Exception as e:
    rep.append(f"ieugwaspy import FAIL: {e}")
try:
    import pandas as pd
    rep.append(f"pandas {pd.__version__}")
except Exception as e:
    rep.append(f"pandas import FAIL: {e}")
open(r'D:\endometriosis_project\_opengwas_check.txt','w',encoding='utf-8').write("\n".join(rep))
print("ok")
