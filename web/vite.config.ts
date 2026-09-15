import { fileURLToPath, URL } from "node:url";
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";

/**
 * A plain single-page application. No server framework, no authentication, no
 * database: the viewer holds no accounts and no session, which is a requirement
 * of the tool rather than a simplification of it — a prediction endpoint has
 * nothing to authenticate, and a public server that makes users register is not
 * eligible for a web-server issue.
 *
 * `/api`, `/evidence`, `/predict`, `/scan` and `/batch` are proxied to the
 * QuadCond service during development so the dev server and a production
 * deployment present the same single origin.
 */
const API = process.env.QUADCOND_URL ?? "http://127.0.0.1:8765";
const proxied = ["/evidence", "/predict", "/scan", "/batch", "/info", "/health", "/ready"];

export default defineConfig({
  plugins: [react(), tailwindcss()],
  resolve: {
    alias: { "@": fileURLToPath(new URL("./src", import.meta.url)) },
  },
  server: {
    port: 8080,
    proxy: Object.fromEntries(
      proxied.map((path) => [path, { target: API, changeOrigin: true }]),
    ),
  },
  build: { outDir: "dist", sourcemap: false },
});
