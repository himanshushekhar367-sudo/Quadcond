# Web-server figures

Three figures for the web-server description of QuadCond / AENNA-3D. Each is in
PNG and LZW-compressed TIFF at 600 dpi, and in PDF with embedded fonts. Widths:
WS1 and WS2 are 180 mm (double column), WS3 is 86 mm (single column).

| Figure | File | What it shows |
| --- | --- | --- |
| WS1 | `Figure_WS1_architecture` | One process serving the workbench and the API; request checks; 12 heads by label type; applicability gate; provenance; verified release assets; optional AlphaGenome join |
| WS2 | `Figure_WS2_interface` | The interface, A–F: layout, K⁺ versus Na⁺ evidence cards, an i-motif card, a refusal, and a condition sweep with a verdict per head |
| WS3 | `Figure_WS3_condition_sweep` | Predicted Tm across 0–150 mM K⁺ or Na⁺ for the telomeric repeat, from `/scan/conditions` |

Legends are in `src/LEGENDS.md`.

## Rebuilding

Everything comes from a running server, so start one first
(`quadcond serve --port 8765`), then from `src/`:

```bash
node capture_screens.js http://127.0.0.1:8765/ raw/   # screenshots for WS2 (needs Playwright)
sh fetch_sweeps.sh http://127.0.0.1:8765              # JSON for WS3
mkdir -p out && python3 ws1.py && python3 ws2.py && python3 ws3.py
```

`ws3.py` refuses to draw if the server's provenance block says the model does not
match the release manifest, or if any swept point is outside the head's domain.
Figure 1 of the manuscript uses the same palette and font, set in
`manuscript/build_manuscript.py`.
