"""Download all TCGA-BRCA metadata needed for Option 1 (no slides).

Outputs (data/metadata/):
  gdc_slides.csv                 every TCGA-BRCA slide on GDC (name, size, DX/TS/BS)
  clinical_patient_brca.txt      GDC BCR Biotab: IHC ER/PR/HER2, HER2 score, FISH, histology, race
  biospecimen_slide_brca.txt     GDC Biotab: % tumour / stroma / necrosis per frozen slide
  pam50.csv                      PAM50 subtype (cBioPortal PanCancer Atlas)
  molecular.csv                  ESR1/PGR/ERBB2 mRNA (RSEM), GISTIC CNA, RPPA (cBioPortal)

Usage:  python scripts/fetch_metadata.py [--out data/metadata]
"""
import argparse
import os

import pandas as pd
import requests

GDC = "https://api.gdc.cancer.gov"
CBIO = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
BIOTAB = {  # GDC file ids (open access)
    "clinical_patient_brca.txt": "8162d394-8b64-4da2-9f5b-d164c54b9608",
    "biospecimen_slide_brca.txt": "c10090b8-fa42-4e32-86a0-9dedceb5b605",
}


def download(file_id, path):
    # GDC sometimes truncates a download; retry until the size matches Content-Length
    for _ in range(3):
        r = requests.get(f"{GDC}/data/{file_id}", timeout=120)
        expected = int(r.headers.get("Content-Length", len(r.content)))
        if len(r.content) == expected:
            with open(path, "wb") as f:
                f.write(r.content)
            return len(r.content)
    raise RuntimeError(f"incomplete download for {file_id}")


def slides(out):
    filt = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-BRCA"}},
        {"op": "=", "content": {"field": "data_type", "value": "Slide Image"}}]}
    r = requests.post(f"{GDC}/files", json={
        "filters": filt, "size": 10000,
        "fields": "file_id,file_name,file_size,experimental_strategy,cases.submitter_id,access"})
    rows = [{"file_id": h["file_id"], "file_name": h["file_name"], "size": h["file_size"],
             "strategy": h.get("experimental_strategy"), "case": h["cases"][0]["submitter_id"],
             "access": h["access"]} for h in r.json()["data"]["hits"]]
    pd.DataFrame(rows).to_csv(os.path.join(out, "gdc_slides.csv"), index=False)
    return len(rows)


def cbioportal(out):
    r = requests.get(f"{CBIO}/studies/{STUDY}/clinical-data",
                     params={"clinicalDataType": "PATIENT", "projection": "SUMMARY", "pageSize": 100000})
    cd = pd.DataFrame(r.json())
    cd[cd.clinicalAttributeId == "SUBTYPE"][["patientId", "value"]].to_csv(
        os.path.join(out, "pam50.csv"), index=False)
    genes = {2099: "ESR1", 5241: "PGR", 2064: "ERBB2"}
    parts = []
    for prof in ["rna_seq_v2_mrna", "gistic", "rppa"]:
        r = requests.post(f"{CBIO}/molecular-profiles/{STUDY}_{prof}/molecular-data/fetch",
                          params={"projection": "SUMMARY"},
                          json={"entrezGeneIds": list(genes), "sampleListId": f"{STUDY}_all"})
        d = pd.DataFrame(r.json())
        d["gene"] = d.entrezGeneId.map(genes)
        d["prof"] = prof
        parts.append(d[["sampleId", "patientId", "gene", "prof", "value"]])
    pd.concat(parts).to_csv(os.path.join(out, "molecular.csv"), index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/metadata")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print("slides:", slides(a.out))
    for name, fid in BIOTAB.items():
        print(name, download(fid, os.path.join(a.out, name)), "bytes")
    cbioportal(a.out)
    print("done ->", a.out)
