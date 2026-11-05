#!/usr/bin/env sh
# Data for Figure WS3: two /scan/conditions calls against a running server.
#   sh fetch_sweeps.sh http://127.0.0.1:8765
BASE=${1:-http://127.0.0.1:8765}
SEQ=GGGTTAGGGTTAGGGTTAGGG
curl -s -X POST "$BASE/scan/conditions" -H 'content-type: application/json' \
  -d "{\"sequence\":\"$SEQ\",\"axis\":\"k\",\"min\":0,\"max\":150,\"points\":16,\"na\":0,\"mg\":0,\"ph\":7.0,\"temperature\":25,\"heads\":[\"g4_tm\"]}" > raw/sweep_k.json
curl -s -X POST "$BASE/scan/conditions" -H 'content-type: application/json' \
  -d "{\"sequence\":\"$SEQ\",\"axis\":\"na\",\"min\":0,\"max\":150,\"points\":16,\"k\":0,\"mg\":0,\"ph\":7.0,\"temperature\":25,\"heads\":[\"g4_tm\"]}" > raw/sweep_na.json
