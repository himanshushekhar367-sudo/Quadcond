# QuadCond / AENNA-3D: one image that serves the HTTP API and the workbench.
#
# The model files are fetched from the GitHub release and checked against
# quadcond/assets_manifest.json while the image is BUILT. A file whose SHA-256
# does not match fails the build, so a running container can only ever hold
# the released model, and needs no network or writable asset directory.
#
#   docker build -t quadcond .
#   docker run -p 7860:7860 quadcond          # then open http://localhost:7860/
#
# The port comes from $PORT (Cloud Run and similar set it); 7860 otherwise.

# ---- 1. build the AENNA-3D workbench --------------------------------------
FROM node:22-slim AS web
WORKDIR /src/web
COPY web/package.json web/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY web/ ./
RUN npm run build

# ---- 2. the service ---------------------------------------------------------
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    QUADCOND_CACHE=/opt/quadcond/assets \
    QUADCOND_WEB_ROOT=/opt/quadcond/web \
    PORT=7860

WORKDIR /opt/quadcond/src
COPY pyproject.toml README.md LICENSE ./
COPY quadcond/ quadcond/
RUN pip install --no-cache-dir .

# Download the released assets and verify every one. `status` exits non-zero
# if a required asset is missing or fails its checksum, which stops the build.
RUN quadcond assets fetch && quadcond assets status

COPY --from=web /src/web/dist /opt/quadcond/web

RUN useradd --create-home --uid 1000 app && chown -R app:app /opt/quadcond
USER app
WORKDIR /home/app

EXPOSE 7860
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import os,urllib.request,sys; sys.exit(0 if urllib.request.urlopen(f'http://127.0.0.1:{os.environ[\"PORT\"]}/ready', timeout=4).status == 200 else 1)"

CMD ["sh", "-c", "exec quadcond serve --host 0.0.0.0 --port \"$PORT\""]
