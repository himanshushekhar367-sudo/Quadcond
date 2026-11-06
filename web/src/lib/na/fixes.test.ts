/**
 * Pins the 0.5.2 workbench fixes: typing into number boxes, refusals that are
 * not outages, a selection that survives a buffer change, and the suggested
 * training buffer of a hidden card.
 */
import { describe, it } from "node:test";
import assert from "node:assert/strict";
import { parseDraft } from "./parse.ts";
import { keepsSelection, REQUEST_DEBOUNCE_MS } from "./store.ts";
import { describeBuffer, isRefusalStatus } from "../quadcond/client.ts";

describe("number boxes accept typing through an out-of-range prefix", () => {
  it("does not read '2' as a temperature when the minimum is 4", () => {
    assert.equal(parseDraft("2", 4, 95), null);
    assert.equal(parseDraft("25", 4, 95), 25);
  });
  it("treats a half-typed decimal and an empty box as no value", () => {
    assert.equal(parseDraft("", 4, 8.5), null);
    assert.equal(parseDraft("6.", 4, 8.5), 6);
    assert.equal(parseDraft("6.5", 4, 8.5), 6.5);
    assert.equal(parseDraft("abc", 0, 10), null);
  });
});

describe("a refused input is not an outage", () => {
  it("classifies 400/413/422 as refusals and 5xx as unavailable", () => {
    for (const s of [400, 413, 422]) assert.equal(isRefusalStatus(s), true);
    for (const s of [500, 502, 503, 404]) assert.equal(isRefusalStatus(s), false);
  });
});

describe("the 3D selection survives a buffer change", () => {
  const members = [{ id: "g-quadruplex-parallel@+" }, { id: "i-motif@-" }];
  it("keeps the selection when only the conditions changed", () => {
    assert.equal(keepsSelection("i-motif@-", "DNA:GGG", "DNA:GGG", members), true);
  });
  it("resets it for a new sequence or a structure that no longer exists", () => {
    assert.equal(keepsSelection("i-motif@-", "DNA:GGG", "DNA:CCC", members), false);
    assert.equal(keepsSelection("duplex@+", "DNA:GGG", "DNA:GGG", members), false);
    assert.equal(keepsSelection("", "DNA:GGG", "DNA:GGG", members), false);
  });
  it("waits briefly before asking, instead of once per pointer move", () => {
    assert.ok(REQUEST_DEBOUNCE_MS >= 100);
    assert.ok(REQUEST_DEBOUNCE_MS <= 400);
  });
});

describe("a hidden card can offer its model's buffer", () => {
  it("describes the suggestion in panel order with units", () => {
    assert.equal(describeBuffer({ temperature: 25, mg: 0, k: 100 }), "K⁺ 100 mM · Mg²⁺ 0 mM · 25 °C");
    assert.equal(describeBuffer({ ph: 7 }), "pH 7");
  });
});
