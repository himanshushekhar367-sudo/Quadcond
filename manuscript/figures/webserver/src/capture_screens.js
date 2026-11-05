// Screenshots for Figure WS2, taken from a running `quadcond serve`.
//   node capture_screens.js http://127.0.0.1:8765/ raw/
// Needs Playwright with Chromium. 1680x1000 viewport at 3x device scale.
const { chromium } = require('playwright');
const [URL = 'http://127.0.0.1:8765/', OUT = 'raw/'] = process.argv.slice(2);
(async () => {
  const b = await chromium.launch({ args: ['--use-gl=swiftshader', '--enable-unsafe-swiftshader'] });
  const page = async (h = 1000) => {
    const p = await (await b.newContext({ viewport: { width: 1680, height: h }, deviceScaleFactor: 3 })).newPage();
    await p.goto(URL, { waitUntil: 'networkidle' }); await p.waitForTimeout(3000); return p;
  };
  const setv = async (p, label, v) => { const el = p.getByLabel(label, { exact: true }); await el.fill(String(v)); await el.press('Enter'); };
  const cards = p => p.locator('aside ul.space-y-3 > li');
  const record = {};

  // A-C: telomeric G-strand, default cytosol-like buffer, then Na+ only
  let p = await page();
  await p.screenshot({ path: OUT + 'A_overview.png' });
  await cards(p).first().screenshot({ path: OUT + 'C_card_K140.png' });
  record.B = (await cards(p).first().innerText()).match(/[\d.]+ °C\s*\[[^\]]+\]/)[0];
  await setv(p, 'Potassium K+ value', 0); await setv(p, 'Sodium Na+ value', 150); await p.waitForTimeout(3000);
  await cards(p).first().screenshot({ path: OUT + 'D_card_Na150.png' });
  record.C = (await cards(p).first().innerText()).match(/[\d.]+ °C\s*\[[^\]]+\]/)[0];
  await p.close();

  // E-D: telomeric C-strand; refusal at the default buffer, then the i-motif card at the reference buffer
  p = await page();
  await p.selectOption('select[aria-label="Choose an example"]', 'telo-c'); await p.waitForTimeout(3000);
  const notice = p.getByText(/hidden because the query falls outside/i).locator('xpath=ancestor::div[contains(@class,"border")][1]').first();
  await notice.screenshot({ path: OUT + 'F_refusal.png' });
  for (const [k, v] of [['Potassium K+ value', 100], ['Sodium Na+ value', 10], ['Magnesium Mg2+ value', 0], ['pH value', 7], ['Temperature value', 25]]) await setv(p, k, v);
  await p.waitForTimeout(3500);
  const n = await cards(p).count();
  for (let i = 0; i < n; i++) {
    const t = await cards(p).nth(i).innerText();
    if (/i-motif/i.test(t)) { await cards(p).nth(i).click(); await p.waitForTimeout(2000); await cards(p).nth(i).screenshot({ path: OUT + 'E_card_imotif.png' }); record.D = t.match(/pH [\d.]+\s*\[[^\]]+\]/)[0]; break; }
  }
  await p.close();

  // F: Compare > Conditions, K+ swept 0-150 mM on the G-strand
  p = await page(1400);
  await p.getByRole('tab', { name: 'Compare' }).first().click(); await p.waitForTimeout(1200);
  await p.getByText('Conditions', { exact: true }).first().click(); await p.waitForTimeout(1200);
  await p.getByRole('button', { name: 'K⁺', exact: true }).click();
  await p.getByRole('button', { name: /Sweep 0–150 mM/ }).click(); await p.waitForTimeout(6000);
  await p.locator('[class*="lg:grid-cols"] > :nth-child(2)').first().screenshot({ path: OUT + 'F_compare_conditions.png' });
  console.log(JSON.stringify(record));
  await b.close();
})();
