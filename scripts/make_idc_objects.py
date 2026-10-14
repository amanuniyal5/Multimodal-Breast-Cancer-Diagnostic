"""Map TCGA-BRCA slides to their copies in the NCI Imaging Data Commons (IDC) and list the DICOM files of each.

Why: from Ada, GDC serves slides at 0.02-0.3 MB/s per connection, so the 429 GB of frozen slides would take weeks.
IDC mirrors the same TCGA slides (converted to DICOM-WSI, readable by OpenSlide >= 4.0) in public cloud buckets;
Google Cloud Storage gave 1.8 MB/s per connection and 10 MB/s with 16 connections on a compute node (14 Oct 2026).

Needs `pip install idc-index` (run on the laptop; the output CSV is small and committed).
Usage: python scripts/make_idc_objects.py [--manifest data/metadata/frozen_manifest.csv]
Output: data/metadata/idc_slides.csv  (slide_barcode, slide_type, series_uuid, object, size): every TS/BS slide in the
        frozen manifest plus every TCGA-BRCA DX slide (for the extraction-equivalence check against TRACE)
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests
from idc_index import IDCClient

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", default="data/metadata/frozen_manifest.csv")
ap.add_argument("--out", default="data/metadata/idc_slides.csv")
a = ap.parse_args()

c = IDCClient()
c.fetch_index("sm_index")
idc = c.sql_query("""
SELECT s.ContainerIdentifier AS slide_barcode, i.series_aws_url, i.series_size_MB
FROM index i JOIN sm_index s ON i.SeriesInstanceUID = s.SeriesInstanceUID
WHERE i.collection_id = 'tcga_brca'
""")
idc["series_uuid"] = idc.series_aws_url.str.extract(r"s3://idc-open-data/([0-9a-f-]+)/")[0]
idc["slide_type"] = idc.slide_barcode.str.extract(r"-(TS|BS|DX|MS)[0-9A-Z]*$")[0]
assert not idc.slide_barcode.duplicated().any(), "one IDC series per slide expected"
print(f"IDC tcga_brca slide series: {len(idc)}; {idc.slide_type.value_counts().to_dict()}")

m = pd.read_csv(a.manifest)
want = idc[idc.slide_barcode.isin(m.slide_barcode) | (idc.slide_type == "DX")]
missing = set(m.slide_barcode) - set(idc.slide_barcode)
print(f"manifest slides found in IDC: {len(m) - len(missing)} of {len(m)}; plus {int((want.slide_type == 'DX').sum())} DX")


def objects(uuid):
    for attempt in range(6):  # the listing API drops connections now and then
        try:
            r = requests.get("https://storage.googleapis.com/storage/v1/b/idc-open-data/o",
                             params={"prefix": f"{uuid}/", "fields": "items(name,size)"}, timeout=60)
            r.raise_for_status()
            return [(uuid, i["name"], int(i["size"])) for i in r.json().get("items", [])]
        except requests.RequestException:
            if attempt == 5:
                raise
            time.sleep(2 * 2 ** attempt)


with ThreadPoolExecutor(6) as ex:
    rows = [x for lst in ex.map(objects, want.series_uuid) for x in lst]
o = pd.DataFrame(rows, columns=["series_uuid", "object", "size"])
out = want[["slide_barcode", "slide_type", "series_uuid"]].merge(o, on="series_uuid")
assert out.series_uuid.nunique() == len(want), "some series have no objects"
out.sort_values(["slide_barcode", "size"], ascending=[True, False]).to_csv(a.out, index=False)
tot = out.groupby("slide_type")["size"].sum() / 1e9
print(f"saved {a.out}: {out.slide_barcode.nunique()} slides, {len(out)} DICOM files; GB by type {tot.round(0).to_dict()}")
