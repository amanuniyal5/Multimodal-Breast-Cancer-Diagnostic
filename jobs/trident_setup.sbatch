#!/bin/bash
# One-time: add TRIDENT (slide segmentation, patching, patch encoders) + OpenSlide to the venv and print its API.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J trident_setup
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:45:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
pip install --no-cache-dir -q openslide-python openslide-bin timm h5py
[ -d $HOME/src/TRIDENT ] || git clone -q https://github.com/mahmoodlab/TRIDENT.git $HOME/src/TRIDENT
cd $HOME/src/TRIDENT && echo "TRIDENT commit $(git rev-parse --short HEAD) $(git log -1 --format=%cd)"
pip install --no-cache-dir -q -e .
cd $HOME/project
python - <<'PY'
import inspect, trident, openslide, timm, torch
print("trident", getattr(trident, "__version__", "?"), "| openslide", openslide.__library_version__, "| timm", timm.__version__, "| torch", torch.__version__)
from trident import OpenSlideWSI
for m in ["segment_tissue", "extract_tissue_coords", "extract_patch_features"]:
    print(m, inspect.signature(getattr(OpenSlideWSI, m)))
from trident.segmentation_models import segmentation_model_factory
print("seg factory", inspect.signature(segmentation_model_factory))
from trident.patch_encoder_models import encoder_factory
print("encoder factory", inspect.signature(encoder_factory))
import trident.patch_encoder_models.load as L
src = inspect.getsource(L)
i = src.find("class UNIInferenceEncoder"); print(src[i:i+2500])
i = src.find("class CustomInferenceEncoder"); print(src[i:i+1500])
PY
echo TRIDENT_OK
