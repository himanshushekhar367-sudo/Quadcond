/**
 * Resolve extensionless relative imports to their .ts / .tsx files.
 *
 * `node --experimental-strip-types` runs TypeScript directly but keeps Node's
 * ESM resolver, which requires a real filename -- while the source uses
 * bundler-style extensionless specifiers (`./geometry`) because Vite resolves
 * them. Rather than rewrite every import in the app to suit the test runner, or
 * pull in a bundler to run three files, this hook fills the gap: try the
 * specifier as given, and on failure try it with each TypeScript extension.
 *
 * Registered by scripts/run-ts-tests.mjs. It affects nothing that ships.
 */
import { existsSync } from "node:fs";
import { fileURLToPath, pathToFileURL } from "node:url";

const CANDIDATES = [".ts", ".tsx", "/index.ts", "/index.tsx"];

export async function resolve(specifier, context, next) {
  try {
    return await next(specifier, context);
  } catch (err) {
    if (!specifier.startsWith(".") && !specifier.startsWith("/")) throw err;
    const base = context.parentURL ? new URL(specifier, context.parentURL) : pathToFileURL(specifier);
    for (const ext of CANDIDATES) {
      const candidate = new URL(base.href + ext);
      if (existsSync(fileURLToPath(candidate))) {
        // No `format` here on purpose. Naming it "module" tells Node the file
        // is already plain JavaScript and skips the type stripper, which then
        // fails on the first `import type`. Letting Node infer from the .ts
        // extension is what routes it through the stripper.
        return { url: candidate.href, shortCircuit: true };
      }
    }
    throw err;
  }
}
