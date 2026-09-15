/**
 * Entry point for `--import`, so the resolver hook is installed in the test
 * runner's child processes too. `node --test` spawns one process per file, and
 * a hook registered only in the parent never reaches them -- which presents as
 * the original "cannot find module" and looks like the hook does not work.
 */
import { register } from "node:module";
register("./ts-extension-resolver.mjs", import.meta.url);
