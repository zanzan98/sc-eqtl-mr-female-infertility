# -*- coding: utf-8 -*-
"""s35b_opentargets_ext.py —— 补充：化学探针 + 通路 + 大疾病列表（生殖/自身免疫筛选）

依赖 s35 已核验的基因 ID。输出：
  tables/41b_opentargets_chemical_probes.csv
  tables/41c_opentargets_pathways.csv
  tables/41d_opentargets_diseases_full.csv （top 100 关联疾病 + 治疗领域）
"""
import json, urllib.request, os, csv

API = "https://api.platform.opentargets.org/api/v4/graphql"
T = r"D:\endometriosis_project\11_sc_eqtl_mr_project\tables"
META = r"D:\endometriosis_project\11_sc_eqtl_mr_project\00_metadata"

GENES = {"CDC42": "ENSG00000070831", "WNT4": "ENSG00000162552", "LINC00339": "ENSG00000218510"}

def gql(query):
    data = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(API, data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

Q = """
{
  target(ensemblId: "%s") {
    approvedSymbol
    chemicalProbes { drugId mechanismOfAction isHighQuality origin probesDrugsScore }
    pathways { pathwayId pathway topLevelTerm }
    associatedDiseases(page: {index: 0, size: 100}) {
      count
      rows { score disease { id name therapeuticAreas { id name } } }
    }
  }
}
"""

probes, paths, diseases = [], [], []
for sym, eid in GENES.items():
    r = gql(Q % eid)
    t = r.get("data", {}).get("target", {})
    if not t:
        print(sym, "NO DATA:", r); continue
    with open(os.path.join(META, f"_ot_ext_{sym}.json"), "w", encoding="utf-8") as f:
        json.dump(r, f, ensure_ascii=False, indent=1)
    for p in (t.get("chemicalProbes") or []):
        probes.append([sym, p.get("drugId"), p.get("mechanismOfAction"),
                       p.get("isHighQuality"), p.get("origin"), p.get("probesDrugsScore")])
    for pw in (t.get("pathways") or []):
        paths.append([sym, pw.get("pathwayId"), pw.get("topLevelTerm"), pw.get("pathway")])
    for row in t.get("associatedDiseases", {}).get("rows", []):
        d = row.get("disease") or {}
        tas = "; ".join(a.get("name","") for a in (d.get("therapeuticAreas") or []))
        diseases.append([sym, d.get("id"), d.get("name"), row.get("score"), tas])
    print(f"{sym}: probes={len(t.get('chemicalProbes') or [])}, pathways={len(t.get('pathways') or [])}, diseases={len(t.get('associatedDiseases',{}).get('rows',[]))}")

def dump(fn, hdr, rows):
    with open(os.path.join(T, fn), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)

dump("41b_opentargets_chemical_probes.csv",
     ["gene","drugId","mechanismOfAction","isHighQuality","origin","probesDrugsScore"], probes)
dump("41c_opentargets_pathways.csv",
     ["gene","pathwayId","topLevelTerm","pathway"], paths)
dump("41d_opentargets_diseases_full.csv",
     ["gene","disease_id","disease_name","overall_score","therapeutic_areas"], diseases)

print(f"\nprobes={len(probes)}, pathways={len(paths)}, diseases={len(diseases)}")

# 筛选生殖/女性健康/自身免疫领域
KW = ["female", "reproduct", "gynecolog", "endometrio", "menstrua", "ovari", "uter",
      "infertil", "pregnan", "testis", "prostate", "breast", "autoimmun", "immune",
      "lupus", "arthritis", "thyroid", "scleroderma", "diabetes"]
print("\n=== 生殖/女性健康/自身免疫相关疾病（关键词命中）===")
for row in diseases:
    txt = (str(row[2]) + " " + str(row[4])).lower()
    if any(k in txt for k in KW):
        print(f"  {row[0]}\t{row[2]}\t score={row[3]:.3g}\t[{row[4]}]")
print("DONE")
