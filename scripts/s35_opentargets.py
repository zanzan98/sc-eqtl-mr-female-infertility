# -*- coding: utf-8 -*-
"""s35_opentargets.py —— 下游：Open Targets 药物靶点 + 简化版 PheWAS（关联疾病谱）

查询 CDC42 / WNT4 / LINC00339 三个基因：
  1) 已知药物与临床候选（drugAndClinicalCandidates）
  2) 关联疾病（associatedDiseases, top 25）+ 治疗领域
  3) 安全风险（safetyLiabilities）
输出 CSV 到 tables/，原始 JSON 到 00_metadata/。

基因 ID（已核验）：
  CDC42 = ENSG00000070831
  WNT4  = ENSG00000162552
  LINC00339 = ENSG00000218510
"""
import json, urllib.request, os, time

API = "https://api.platform.opentargets.org/api/v4/graphql"
T = r"D:\endometriosis_project\11_sc_eqtl_mr_project\tables"
META = r"D:\endometriosis_project\11_sc_eqtl_mr_project\00_metadata"
os.makedirs(META, exist_ok=True)

GENES = {
    "CDC42":    "ENSG00000070831",
    "WNT4":     "ENSG00000162552",
    "LINC00339":"ENSG00000218510",
}

def gql(query):
    data = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(API, data=data,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

Q_TMPL = """
{
  target(ensemblId: "%s") {
    id
    approvedSymbol
    approvedName
    biotype
    functionDescriptions
    tractability { label modality value }
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug { id name drugType maximumClinicalStage }
        diseases { diseaseFromSource }
      }
    }
    associatedDiseases(page: {index: 0, size: 25}) {
      count
      rows {
        score
        disease { id name therapeuticAreas { id name } }
        datatypeScores { id score }
      }
    }
    safetyLiabilities { event effects { direction dosing } datasource url }
  }
}
"""

results = {}
for sym, eid in GENES.items():
    q = Q_TMPL % eid
    try:
        r = gql(q)
        results[sym] = r
        with open(os.path.join(META, f"_ot_{sym}.json"), "w", encoding="utf-8") as f:
            json.dump(r, f, ensure_ascii=False, indent=1)
        t = r.get("data", {}).get("target")
        nd = t.get("drugAndClinicalCandidates", {}).get("count") if t else None
        na = t.get("associatedDiseases", {}).get("count") if t else None
        print(f"{sym}: drugs={nd}, assocDiseases={na}")
    except Exception as e:
        print(f"{sym}: ERROR {e}")
        results[sym] = {"error": str(e)}
    time.sleep(1)

# ---- 汇总药物表 ----
import csv
drug_rows = []
for sym, r in results.items():
    t = r.get("data", {}).get("target")
    if not t: continue
    for row in t.get("drugAndClinicalCandidates", {}).get("rows", []):
        d = row.get("drug") or {}
        dis = row.get("diseases") or []
        dis_names = "; ".join(x.get("diseaseFromSource","") for x in dis if isinstance(x, dict))
        drug_rows.append([
            sym, d.get("id"), d.get("name"), d.get("drugType"),
            d.get("maximumClinicalStage"), row.get("maxClinicalStage"), dis_names
        ])
with open(os.path.join(T, "41_opentargets_drugs.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene","drug_id","drug_name","drug_type","drug_max_clinical_stage",
                "target_max_stage","diseases"])
    w.writerows(drug_rows)

# ---- 汇总疾病表 ----
dis_rows = []
for sym, r in results.items():
    t = r.get("data", {}).get("target")
    if not t: continue
    for row in t.get("associatedDiseases", {}).get("rows", []):
        dis = row.get("disease") or {}
        tas = "; ".join(a.get("name","") for a in (dis.get("therapeuticAreas") or []))
        ds = "; ".join(f"{x.get('id')}={x.get('score'):.3g}" for x in (row.get("datatypeScores") or []) if x.get("score") is not None)
        dis_rows.append([sym, dis.get("id"), dis.get("name"), row.get("score"), tas, ds])
with open(os.path.join(T, "41_opentargets_diseases.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene","disease_id","disease_name","overall_score","therapeutic_areas","datatype_scores"])
    w.writerows(dis_rows)

# ---- 安全风险 ----
safe_rows = []
for sym, r in results.items():
    t = r.get("data", {}).get("target")
    if not t: continue
    for s in (t.get("safetyLiabilities") or []):
        eff = s.get("effects") or []
        eff_s = "; ".join(f"{x.get('direction','')}/{x.get('dosing','')}" for x in eff if isinstance(x, dict))
        safe_rows.append([sym, s.get("event"), eff_s, s.get("datasource"), s.get("url")])
with open(os.path.join(T, "41_opentargets_safety.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene","event","effects(direction/dosing)","datasource","url"])
    w.writerows(safe_rows)

print(f"\n药物行={len(drug_rows)}, 疾病行={len(dis_rows)}, 安全行={len(safe_rows)}")
print("DONE")
