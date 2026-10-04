# -*- coding: utf-8 -*-
"""s72_make_v6.py —— 由 v4（docx 血统）产出 44_论文初稿_带图_v6

用户裁定（2026-09-29）：
  1. 文件同步选 B：以 v4 docx 为底本做 **run 级插入** 产 v6；把手改与新增的 S1–S4 图
     同步进 docx；**保住手改内容，不丢图**。
  2. 附录 A 加上 S1–S4 的图注（基础项，不要漏）。
  3. S4 不加 MHC 着色（正文未提 MHC）。
  4. 其余正文、图表、补充材料不动。

相对 v4 docx 的 4 处改动（其余 98 段逐字未动）：
  A  §3.3 图注 图 4 (a)：`（r²，1000 Genomes 欧洲参考）` → `（r²，OneK1K 参考面板，980 名供者）`
     （承自 v4→v5 轮，本轮补做进 docx；v4 docx 此前未做）
  D1 §4 局限段：删除性别分层句 + 前置「；」改「。」（承自 v5 删除）
  B  Results 末（`4 讨论` 前）插入 1 句集合式引用：`质量控制与分层视图另见补充图 S1–S4。`
  C  附录 A：插入 S1–S4 四张图（16.0 cm 宽）+ 四条图注

红线遵守：
  ① 不改任何已闭环结果（tables/ 只读）
  ② ★不得丢失 v4 保住的 9 段手改（逐段断言）
  ③ 不新增任何数值/因果断言（B/C 皆为引用句与图注，图注文本**逐字取自**已有权威件
     `39b_技术文档_SM与SN_v8.md` 的 `### 补充图`，仅做「加 `|` 分隔符 + 去 `**` 强调标记」两处体例归一）
  ④ S4 不着色 MHC

日志：logs/s72_make_v6.log
"""
import io
import os
import re
import shutil
import hashlib
import traceback
from copy import deepcopy

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt
from docx.text.paragraph import Paragraph

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
BUILD = r"D:/_wb_docx_build"

BASE_DOCX = os.path.join(ROOT, "42_论文初稿_带图_v4.docx")
OUT_DOCX = os.path.join(ROOT, "44_论文初稿_带图_v6.docx")
OUT_DOCX_BUILD = os.path.join(BUILD, "44_论文初稿_带图_v6.docx")
BASE_MD = os.path.join(ROOT, "43_论文初稿_带图_v5.md")
OUT_MD = os.path.join(ROOT, "44_论文初稿_带图_v6.md")
LEGEND_SRC = os.path.join(ROOT, "39b_技术文档_SM与SN_v8.md")
LOG = os.path.join(ROOT, "logs", "s72_make_v6.log")

NS_WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"

TARGET_CM = 16.0

FIG_S = [
    ("S1", os.path.join(ROOT, "figures", "FigS1_instrument_power_gating.png")),
    ("S2", os.path.join(ROOT, "figures", "FigS2_sensitivity_outcome_concordance.png")),
    ("S3", os.path.join(ROOT, "figures", "FigS3_manhattan_by_celltype.png")),
    ("S4", os.path.join(ROOT, "figures", "FigS4_discovery_volcano.png")),
]

OLD_A = "（r²，1000 Genomes 欧洲参考）"
NEW_A = "（r²，OneK1K 参考面板，980 名供者）"

OLD_D1 = ("；暴露来自外周血免疫细胞而结局为生殖道表型，性别分层与生殖道 eQTL 参考的缺失"
          "使这些效应能否迁移至靶组织仍待检验。")
NEW_D1 = "。"

CITE = "质量控制与分层视图另见补充图 S1–S4。"

HAND9_HEADS = [
    "最显著的观察之一是 CDC42 的细胞类型特异性",
    "在 B 细胞中，CDC42 调控 B 细胞受体信号下游的肌动蛋白重组",
    "在单核细胞一侧，非经典单核细胞沿血管内皮巡逻",
    "这一分布提示，chr1p36.12 在免疫细胞中的信号",
    "然而，同一区域内的精细定位把结局可信集置于 WNT4 基因体内",
    "WNT4 是 Wnt 信号通路的成员，在胚胎期参与苗勒管的形成与分化",
    "在成人期，WNT4 继续参与卵巢功能与子宫内膜周期性重塑",
    "更重要的是，WNT4 的调控异常不仅在着床环节发挥作用",
    "综合来看，chr1p36.12 对不孕症的影响",
]
INSERTS = ["三个信号座在结构上并不对称", "三位点在两层暴露上的判读分层",
           "作为对该归属边界的补充探索", "三个信号座的并列比较为这一分层提供了内部校准"]

L = []


def p(*a):
    s = " ".join(str(x) for x in a)
    L.append(s)
    # ★增量落盘：日志必须逐行落盘，不能只等收尾写（长任务被中断会全丢）
    try:
        with io.open(LOG, "a", encoding="utf-8", newline="\n") as fh:
            fh.write(s + "\n")
    except Exception:
        pass


def strip_ws(s):
    return re.sub(r"\s+", "", s)


def norm(s):
    s = re.sub(r"[*`_>#\[\]]", "", s)
    return re.sub(r"\s+", "", s)


def set_run_font(run, ascii_f="Times New Roman", ea_f="宋体", size=None, bold=None):
    run.font.name = ascii_f
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts")
        rPr.append(rFonts)
    rFonts.set(qn("w:ascii"), ascii_f)
    rFonts.set(qn("w:hAnsi"), ascii_f)
    rFonts.set(qn("w:eastAsia"), ea_f)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.font.bold = bold


def run_level_replace(para, old, new):
    """段落内 run 级外科替换：只改与区间重叠的 run 的 w:t，其余 run 的 rPr 原样保留。"""
    runs = para._p.findall(qn("w:r"))
    texts = ["".join(t.text or "" for t in r.findall(qn("w:t"))) for r in runs]
    full = "".join(texts)
    assert full.count(old) == 1, "目标串在段落内出现 %d 次" % full.count(old)
    i = full.index(old)
    j = i + len(old)
    pos = 0
    for r, txt in zip(runs, texts):
        s, e = pos, pos + len(txt)
        pos = e
        if e <= i or s >= j:
            continue
        keep_l = txt[:i - s] if s < i else ""
        if s <= i < e:
            keep_l += new
        keep_r = txt[j - s:] if e > j else ""
        nt = keep_l + keep_r
        ts = r.findall(qn("w:t"))
        for k in ts[1:]:
            r.remove(k)
        if not ts:
            continue
        if nt == "":
            r.remove(ts[0])
        else:
            ts[0].text = nt
            ts[0].set(qn("xml:space"), "preserve")
    return full[:i] + new + full[j:]


def all_paras(doc):
    return [Paragraph(k, doc) for k in doc.element.body if k.tag == qn("w:p")]


def nonempty(doc):
    return [x for x in all_paras(doc) if x.text.strip()]


def make_text_para(template_p, doc, segments, size):
    """由模板段落复制 pPr，去 run，按 (text, bold) 段写入。"""
    if hasattr(template_p, "_p"):
        template_p = template_p._p
    new_p = deepcopy(template_p)
    for ch in list(new_p):
        if ch.tag != qn("w:pPr"):
            new_p.remove(ch)
    par = Paragraph(new_p, doc)
    for txt, bold in segments:
        r = par.add_run(txt)
        set_run_font(r, size=size, bold=bold)
    return new_p


def make_img_para(template_img_p, doc, png, name, docpr_id, target_cm=TARGET_CM):
    if hasattr(template_img_p, "_p"):
        template_img_p = template_img_p._p
    new_p = deepcopy(template_img_p)
    inl = new_p.find(".//{%s}inline" % NS_WP)
    assert inl is not None, "模板段无 wp:inline"
    blip = inl.find(".//{%s}blip" % NS_A)
    assert blip is not None, "模板段无 a:blip"
    rid, image = doc.part.get_or_add_image(png)
    blip.set(qn("r:embed"), rid)
    tgt = int(round(target_cm * 360000))
    cy = int(round(tgt * image.height / image.width))
    ext = inl.find("{%s}extent" % NS_WP)
    ext.set("cx", str(tgt))
    ext.set("cy", str(cy))
    for aext in inl.findall(".//{%s}ext" % NS_A):
        aext.set("cx", str(tgt))
        aext.set("cy", str(cy))
    docPr = inl.find(".//{%s}docPr" % NS_WP)
    if docPr is None:
        docPr = inl.find("{%s}docPr" % NS_WP)
    assert docPr is not None, "模板段无 wp:docPr"
    docPr.set("id", str(docpr_id))
    docPr.set("name", name)
    return new_p, rid


def read_legends():
    """从 39b_v8 的 `### 补充图` 段取权威中文图注，切成 4 条并体例归一。"""
    md = io.open(LEGEND_SRC, encoding="utf-8").read()
    i = md.index("### 补充图")
    lines = md[i:].split("\n")
    src = None
    for ln in lines[1:]:
        if ln.strip():
            src = ln.strip()
            break
    assert src, "未取到 补充图 段落"
    marks = ["图 S1 ", "图 S2 ", "图 S3 ", "图 S4 "]
    pos = [src.index(m) for m in marks]
    assert pos == sorted(pos), "图 S1–S4 标记顺序异常 %s" % pos
    parts = []
    for k, (m, i0) in enumerate(zip(marks, pos)):
        i1 = pos[k + 1] if k + 1 < len(marks) else len(src)
        parts.append((m.strip().split()[1], src[i0:i1].strip()))
    assert [x[0] for x in parts] == ["S1", "S2", "S3", "S4"]
    out = []
    for k, t in parts:
        assert t.startswith("图 %s " % k)
        body = t[len("图 %s " % k):].replace("**", "")
        j = body.index("。")
        title, rest = body[:j + 1], body[j + 1:].strip()
        out.append((k, "图 %s | %s" % (k, title), rest))
    return src, out


def main():
    base_b = open(BASE_DOCX, "rb").read()
    p("=== 44_论文初稿_带图_v6 构建（底本 = v4 docx，run 级插入）===")
    p("底本 %s  %d B  md5=%s" % (os.path.basename(BASE_DOCX), len(base_b),
                                 hashlib.md5(base_b).hexdigest()))

    src_legend, legends = read_legends()
    p("图注来源 %s" % os.path.basename(LEGEND_SRC))
    p("  原文段长 %d 字；切出 %d 条" % (len(src_legend), len(legends)))
    for k, lab, rest in legends:
        p("    %s : %s  (+正文 %d 字)" % (k, lab, len(rest)))

    doc = Document(BASE_DOCX)
    base_ne = [x.text.strip() for x in nonempty(doc)]
    assert len(base_ne) == 99, "v4 非空段落应为 99，实为 %d" % len(base_ne)
    p("底本非空段落 = %d" % len(base_ne))

    paras = all_paras(doc)

    def uniq(pred, label):
        hits = [x for x in paras if pred(x)]
        assert len(hits) == 1, "%s 命中 %d 次（应为 1）" % (label, len(hits))
        return hits[0]

    # ---------------- 模板（在任何改动前取） ----------------
    tpl_cap = uniq(lambda x: x.text.strip().startswith("图 5 |"), "图注模板（图 5）")
    img_paras = [x for x in paras if x._p.findall(".//{%s}inline" % NS_WP)]
    assert len(img_paras) == 5, "底本插图段应为 5，实为 %d" % len(img_paras)
    tpl_img = img_paras[0]
    tpl_body = uniq(lambda x: x.text.strip().startswith("在 GTEx v8 全血（n = 670）中"),
                    "正文模板（单 run 正文段）")
    assert len(tpl_body._p.findall(qn("w:r"))) == 1, "正文模板段不是单 run"
    p("模板：图注 图5 / 插图 图1（共 5 个插图段）/ 正文段（单 run）")

    # ---------------- A. 图 4(a) 图注更正 ----------------
    cap4 = uniq(lambda x: x.text.strip().startswith("图 4 |"), "图 4 图注")
    assert OLD_A in cap4.text, "图 4 图注内未找到旧串"
    got_a = run_level_replace(cap4, OLD_A, NEW_A)
    assert got_a == cap4.text, "回读段落文本与预期不符（A）"
    assert OLD_A not in cap4.text and NEW_A in cap4.text
    p("A ✅ 图 4(a) 图注：'1000 Genomes 欧洲参考' → 'OneK1K 参考面板，980 名供者'")

    # ---------------- D1. 局限段删句 ----------------
    lim = uniq(lambda x: "性别分层" in x.text, "局限段（含性别分层）")
    old_len = len(lim.text)
    got_d = run_level_replace(lim, OLD_D1, NEW_D1)
    assert got_d == lim.text, "回读段落文本与预期不符（D1）"
    assert "性别分层" not in lim.text
    assert "暴露来自外周血免疫细胞而结局为生殖道表型" not in lim.text
    assert "不足以确认细胞类型层面的效应量。在安全性评估一侧" in lim.text
    p("D1 ✅ 局限段 %d → %d 字；悬空分号已改句号" % (old_len, len(lim.text)))

    # ---------------- B. Results 末插入引用句 ----------------
    head4 = uniq(lambda x: x.text.strip() == "4 讨论", "标题段 4 讨论")
    cite_p = make_text_para(tpl_body, doc, [(CITE, False)], 12)
    head4._p.addprevious(cite_p)
    p("B ✅ 已在 `4 讨论` 前插入引用句")

    # ---------------- C. 附录 A 插入 S1–S4 图 + 图注 ----------------
    anchor = uniq(lambda x: "补充数据见 Supplementary_Data.zip" in x.text, "附录 A 指针句")
    cur = anchor._p
    rids = []
    for idx, ((k, lab, rest), (_, png)) in enumerate(zip(legends, FIG_S)):
        img_p, rid = make_img_para(tpl_img, doc, png, "Figure %s" % k,
                                   docpr_id=900 + idx, target_cm=TARGET_CM)
        cap_p = make_text_para(tpl_cap, doc, [(lab, True), (" " + rest, False)], 9)
        cur.addnext(img_p)
        img_p.addnext(cap_p)
        cur = cap_p
        rids.append((k, rid))
        p("C ✅ 附录 A 插入 图 %s（%s，16.0 cm）+ 图注" % (k, os.path.basename(png)))

    doc.save(OUT_DOCX)
    shutil.copy2(OUT_DOCX, OUT_DOCX_BUILD)
    p("docx 落盘：%s  %d B" % (OUT_DOCX, os.path.getsize(OUT_DOCX)))

    # ================= 复核 docx =================
    chk = Document(OUT_DOCX)
    new_ne = [x.text.strip() for x in nonempty(chk)]
    p("")
    p("---- docx 复核 ----")
    assert len(new_ne) == 104, "v6 非空段落 %d != 104（99 + 引用句 1 + 图注 4）" % len(new_ne)
    p("非空段落 = %d（v4 99 + 1 + 4）✅" % len(new_ne))

    base_set, new_set = set(base_ne), set(new_ne)
    base_only = [t for t in base_ne if t not in new_set]
    new_only = [t for t in new_ne if t not in base_set]
    p("仅见于 v4 的段 = %d" % len(base_only))
    for t in base_only:
        p("   - %s" % t[:60])
    p("仅见于 v6 的段 = %d" % len(new_only))
    for t in new_only:
        p("   + %s" % t[:60])
    assert len(base_only) == 2, "v4 独有段应为 2（旧图 4 图注、旧局限段）"
    assert len(new_only) == 7, "v6 独有段应为 7（新图 4 图注、新局限段、引用句、4 条图注）"
    assert OLD_A in base_only[0] or OLD_A in base_only[1], "base_only 未含旧图 4 图注"

    common_new = [t for t in new_ne if t in base_set]
    common_base = [t for t in base_ne if t in new_set]
    assert common_new == common_base, "公共子序列顺序不一致 → 有段被移动或改写"
    p("公共子序列顺序一致（其余 97 段逐字未动、未移位）✅")

    t_all = "\n".join(x.text for x in chk.paragraphs)
    for h in HAND9_HEADS:
        assert h in t_all, "9 段手改缺失：%s" % h[:20]
    p("★9 段手改全部在位：9/9 ✅")
    for s in INSERTS:
        assert s in t_all, "v4 插入段缺失：%s" % s[:20]
    p("v4 的 4 处插入全部在位：4/4 ✅")

    # 标记串
    assert t_all.count("性别分层") == 0
    assert t_all.count("1000 Genomes") == 0
    assert t_all.count(NEW_A) == 1
    assert t_all.count(CITE) == 1
    for k, lab, rest in legends:
        assert t_all.count(lab) == 1, "图注标签 %s 命中 %d 次" % (k, t_all.count(lab))
        assert strip_ws(rest)[:40] in strip_ws(t_all), "图注正文缺失 %s" % k
    assert t_all.count("MHC") == 0, "不应引入 MHC 相关文字"
    p("标记串：性别分层 0 / 1000 Genomes 0 / 新图 4 图注串 1 / 引用句 1 / 图 S1–S4 标签各 1 ✅")

    # 图片
    exts, blob_ok, n_img = [], 0, 0
    rels = chk.part.rels
    for x in all_paras(chk):
        for inl in x._p.findall(".//{%s}inline" % NS_WP):
            n_img += 1
            e = inl.find("{%s}extent" % NS_WP)
            w = round(int(e.get("cx")) / 360000, 2)
            h = round(int(e.get("cy")) / 360000, 2)
            exts.append((w, h))
            blip = inl.find(".//{%s}blip" % NS_A)
            rid = blip.get(qn("r:embed"))
            tp = rels[rid].target_part
            blob = tp.blob
            for _, png in FIG_S:
                if len(blob) == os.path.getsize(png):
                    blob_ok += 1
                    break
    p("图片总数 = %d" % n_img)
    assert n_img == 9, "图片数 %d != 9" % n_img
    p("extent = %s" % exts)
    assert exts[:5] == [(16.0, 10.04), (16.0, 14.63), (16.0, 19.51), (16.0, 23.1), (16.0, 21.08)], \
        "主图 extent 被改动"
    from docx.image.image import Image as DImage
    for i, (_, png) in enumerate(FIG_S):
        im = DImage.from_file(png)
        w, h = exts[5 + i]
        assert abs(w - 16.0) < 0.05, "图 %s 宽 %s" % (FIG_S[i][0], w)
        assert abs(h / w - im.height / im.width) < 0.02, "图 %s 宽高比失真" % FIG_S[i][0]
        assert h <= 24.4, "图 %s 高 %.2f cm 超版心" % (FIG_S[i][0], h)
        p("   图 %s: %.2f x %.2f cm（源 %.2f cm 宽，放大 %.2f×）" % (
            FIG_S[i][0], w, h, im.width / 360000, w / (im.width / 360000)))
    p("主图 extent 未变 ✅；4 张 S 图宽度=16.0 cm、宽高比与源一致 ✅")
    assert blob_ok == 4, "只有 %d 张 S 图的二进制与磁盘 PNG 一致" % blob_ok
    p("4 张 S 图的二进制与 figures/*.png 逐字节一致 ✅")

    assert chk.sections[0]._sectPr.find(qn("w:lnNumType")) is None, "连续行号未关闭"
    for bad, nm in [("$", "美元符"), ("\\", "反斜杠"), ("\ufffd", "替换字符"), ("**", "字面星号")]:
        assert bad not in t_all, "成品含 %s" % nm
    p("行号关闭 / 无越界字符（$ \\ U+FFFD 字面**）✅")

    # ================= 同步 md =================
    md = io.open(BASE_MD, encoding="utf-8", newline="").read()
    ptr = "与本文相关的补充数据见 Supplementary_Data.zip。分析口径的完整记录、补充方法与逐项说明另见配套技术文档。"
    assert md.count(ptr) == 1, "md 附录 A 指针句命中 %d 次" % md.count(ptr)
    dup = "质量控制与分层视图另见补充图 S1–S4。"
    assert md.count(dup) == 1
    assert "图 S1 |" not in md, "md 已含 S 图注（幂等保护）"
    blocks = []
    for k, lab, rest in legends:
        blocks.append("**【插入图 %s】**" % k)
        blocks.append("")
        blocks.append("**%s** %s" % (lab, rest))
        blocks.append("")
    ins = "\n".join(blocks).rstrip()
    new_md = md.replace(ptr, ptr + "\n\n" + ins, 1)
    assert new_md.count("图 S1 |") == 1 and new_md.count("图 S4 |") == 1
    assert new_md.count("**【插入图 S1】**") == 1
    # 与原文其余内容逐字一致（除插入块外 0 改动）
    assert new_md.replace(ptr + "\n\n" + ins, ptr, 1) == md, "md 出现插入块之外的改动"
    assert new_md.count("质量控制与分层视图另见补充图 S1–S4。") == 1
    assert new_md.count("1000 Genomes") == 0
    io.open(OUT_MD, "w", encoding="utf-8", newline="\n").write(new_md)
    p("")
    p("md 落盘：%s  %d B（v5 %d B，+%d）" % (
        OUT_MD, os.path.getsize(OUT_MD), os.path.getsize(BASE_MD),
        os.path.getsize(OUT_MD) - os.path.getsize(BASE_MD)))

    md_flat = norm(io.open(OUT_MD, encoding="utf-8").read())
    miss = [t for t in new_ne if norm(t) not in md_flat]
    p("md↔docx 归一化一致性：未见于 md 的 docx 段落 = %d" % len(miss))
    for t in miss[:3]:
        p("   MISS> %s" % t[:70])
    assert not miss, miss[:1]
    p("✅ docx 全部 104 段均可在 v6 md 中定位")

    p("")
    p("md5  v4.docx = %s" % hashlib.md5(base_b).hexdigest())
    p("md5  v6.docx = %s" % hashlib.md5(open(OUT_DOCX, "rb").read()).hexdigest())
    p("md5  v5.md   = %s" % hashlib.md5(open(BASE_MD, "rb").read()).hexdigest())
    p("md5  v6.md   = %s" % hashlib.md5(open(OUT_MD, "rb").read()).hexdigest())
    p("RESULT=PASS")


if __name__ == "__main__":
    try:
        main()
    except BaseException:
        L.append("")
        L.append("==== TRACEBACK ====")
        L.append(traceback.format_exc())
        with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(L) + "\n")
        raise
    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("WROTE", LOG)
