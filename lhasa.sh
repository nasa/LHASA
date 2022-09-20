conda activate lhasa

python lhasa.py --small -icd 2 -scd 2 -t 4

python pfdf/scripts/request.py --filepath pfdf \
  --firms_path pfdf/firms
python pfdf/scripts/gee_export.py --filepath pfdf --gee_username username \
  --firms_path pfdf/firms
python pfdf/scripts/query.py --filepath pfdf --gee_username username \
  --asset_path pfdf
python pfdf/scripts/pull_imerg.py --path pfdf \
  --output_path pfdf/predictions
