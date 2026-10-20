import { test, expect, Page } from '@playwright/test';
import fs from 'node:fs';
import path from 'node:path';

const root = path.resolve(__dirname, '..');
const EXAMPLE_FILE = path.join(root, 'web', 'examples', 'rig-i-like-receptor-innatedb.xgmml');
const exampleXml = fs.readFileSync(EXAMPLE_FILE, 'utf8');

// Expected numbers, computed from the file itself (not from the app):
const EXAMPLE_NODES = (exampleXml.match(/<node /g) || []).length;
const edgeTags = exampleXml.match(/<edge [^>]*>/g) || [];
const EXAMPLE_EDGES = edgeTags.filter((e) => {
  const s = /source="([^"]*)"/.exec(e)![1];
  const t = /target="([^"]*)"/.exec(e)![1];
  return s !== t; // the parser removes self-edges
}).length;
const MAPK_NODES = (exampleXml.match(/<node label="[^"]*"/g) || []).filter((n) => /mapk/i.test(n)).length;

/** Collects everything that must stay empty: console errors/warnings, uncaught errors, dialogs, foreign origins, HTTP errors. */
function watch(page: Page) {
  const problems: string[] = [];
  const origins = new Set<string>();
  page.on('console', (m) => {
    if (m.type() === 'error' || m.type() === 'warning') problems.push(`console.${m.type()}: ${m.text()}`);
  });
  page.on('pageerror', (e) => problems.push(`pageerror: ${e.message}`));
  page.on('dialog', async (d) => {
    problems.push(`dialog(${d.type()}): ${d.message()}`);
    await d.dismiss();
  });
  page.on('request', (r) => {
    const u = r.url();
    if (u.startsWith('data:') || u.startsWith('blob:')) return;
    origins.add(new URL(u).origin);
  });
  page.on('response', (r) => {
    if (r.status() >= 400) problems.push(`HTTP ${r.status()}: ${r.url()}`);
  });
  return { problems, origins };
}

async function openPanel(page: Page, id: string) {
  if (!(await page.locator(`#${id}.in`).count())) {
    await page.click(`a[data-target="#${id}"]`);
  }
  await page.waitForSelector(`#${id}.in:not(.collapsing)`);
}

async function loadExample(page: Page) {
  await page.click('#load_example');
  await expect(page.locator('#Nnodes')).toHaveText(`${EXAMPLE_NODES} nodes`);
}

const positions = (page: Page) =>
  page.evaluate(() =>
    (window as any).d3.selectAll('.node').data().slice(0, 40).map((d: any) => [d.x, d.y]),
  );

const moved = (a: number[][], b: number[][]) => a.some((p, i) => Math.abs(p[0] - b[i][0]) > 1e-6 || Math.abs(p[1] - b[i][1]) > 1e-6);

test('all pages: no console errors, no foreign origins, no analytics, banner is dismissible', async ({ page }) => {
  const w = watch(page);
  for (const p of ['index.html', 'help.html', 'citation.html', 'contact.html']) {
    await page.goto('/' + p);
    await page.waitForLoadState('networkidle');
    await expect(page.locator('#legacy-banner')).toContainText('CellNetVis 1 (2017)');
    await expect(page.locator('#legacy-banner')).toContainText('https://github.com/heberleh/cellnetvis');
    const html = await page.content();
    expect(html).not.toMatch(/googletagmanager|google-analytics|gtag\(|maxcdn|fonts\.googleapis/);
  }
  // dismiss (for the visit only: nothing is stored)
  await page.goto('/help.html');
  await page.click('#legacy-banner button.close');
  await expect(page.locator('#legacy-banner')).toHaveCount(0);
  expect(await page.evaluate(() => localStorage.length + sessionStorage.length + document.cookie.length)).toBe(0);
  await page.reload();
  await expect(page.locator('#legacy-banner')).toHaveCount(1);

  expect([...w.origins]).toEqual(['http://localhost:' + (process.env.PORT || 8080)]);
  expect(w.problems).toEqual([]);
});

test('convertToCellNetVisCC maps values to cellular components', async ({ page }) => {
  await page.goto('/index.html');
  await expect(page.locator('.group_component')).toHaveCount(20);
  await page.waitForFunction(() => Object.keys((window as any).cell).length >= 20);
  const cases: Record<string, string> = {
    mitochondrial_inner_membrane: 'mitochondrion',
    cell_surface: 'extracellular',
    golgi_membrane: 'golgi_apparatus',
    foo: 'not_found',
    '': 'not_found',
    nucleolus: 'nucleus',
  };
  for (const [input, expected] of Object.entries(cases)) {
    const got = await page.evaluate((v) => (window as any).convertToCellNetVisCC(v), input);
    expect(got, `convertToCellNetVisCC(${JSON.stringify(input)})`).toBe(expected);
  }
  // capitalize is defined once and capitalises every word
  expect(await page.evaluate(() => 'plasma membrane'.replace(/_/g, ' ').capitalize())).toBe('Plasma Membrane');
});

test('smoke flow on the example network', async ({ page }) => {
  const w = watch(page);
  await page.goto('/index.html');

  // the cell diagram
  await expect(page.locator('#diagram')).toBeVisible();
  await expect(page.locator('#diagram .group_component')).toHaveCount(20);

  // load example: node and edge counts
  await loadExample(page);
  await expect(page.locator('#Nedges')).toHaveText(`${EXAMPLE_EDGES} links`);
  console.log(`example: ${EXAMPLE_NODES} nodes, ${EXAMPLE_EDGES} edges (file has ${edgeTags.length} edge tags, self-edges removed)`);
  await expect(page.locator('#diagram .node')).toHaveCount(EXAMPLE_NODES);
  await expect(page.locator('#diagram .edge')).toHaveCount(EXAMPLE_EDGES);
  await expect(page.locator('#loading-container')).toBeHidden();

  // the force layout is running: positions change between two samples 500 ms apart
  const a = await positions(page);
  await page.waitForTimeout(500);
  const b = await positions(page);
  expect(moved(a, b)).toBe(true);

  // Stop, then Start
  await openPanel(page, 'collapseForces');
  await page.click('#stop_forces');
  await page.waitForTimeout(100);
  const c = await positions(page);
  await page.waitForTimeout(500);
  const d = await positions(page);
  expect(moved(c, d)).toBe(false);
  await page.click('#restart_forces');
  const e = await positions(page);
  await page.waitForTimeout(500);
  const f = await positions(page);
  expect(moved(e, f)).toBe(true);

  // network measures
  await openPanel(page, 'collapseNodes');
  await page.click('#netmeas');
  await expect(page.locator('#color option[value="[degree]"]')).toHaveCount(1);
  await expect(page.locator('#change_label option[value="[degree]"]')).toHaveCount(1);
  await page.click('#netmeas'); // computing again must not duplicate the options
  await page.waitForTimeout(500);
  await expect(page.locator('#color option[value="[degree]"]')).toHaveCount(1);
  await expect(page.locator('#color option[value="[betweenness]"]')).toHaveCount(1);

  // labels
  const visibleLabels = () =>
    page.evaluate(() => Array.from(document.querySelectorAll('.node_label')).filter((l) => (l as SVGElement).style.visibility === 'visible').length);
  expect(await visibleLabels()).toBe(0);
  await page.click('#showlabels');
  expect(await visibleLabels()).toBe(EXAMPLE_NODES);
  await page.click('#hidelabels');
  expect(await visibleLabels()).toBe(0);
  await page.click('#showlabels');
  expect(await visibleLabels()).toBe(EXAMPLE_NODES);
  await openPanel(page, 'collapseSearch');
  await page.click('#clear'); // hides the labels again
  expect(await visibleLabels()).toBe(0);
  await expect(page.locator('#color')).toHaveValue('Selected CC');

  // search "MAPK", Filter, Clear
  expect(MAPK_NODES).toBeGreaterThan(0);
  await openPanel(page, 'collapseSearch');
  await page.fill('#search_text', 'MAPK');
  const selected = () => page.evaluate(() => (window as any).d3.selectAll('.node').data().filter((d: any) => d.selected).length);
  expect(await selected()).toBe(MAPK_NODES);
  await page.click('#keep');
  const hidden = () =>
    page.evaluate(() => Array.from(document.querySelectorAll('.node')).filter((n) => (n as SVGElement).style.visibility === 'hidden').length);
  expect(await hidden()).toBe(EXAMPLE_NODES - MAPK_NODES);
  await page.click('#clear');
  expect(await hidden()).toBe(0);
  expect(await selected()).toBe(0);
  await page.fill('#search_text', '');

  // donut chart
  await openPanel(page, 'collapseDistribution');
  await expect.poll(() => page.locator('#donutchart path').count()).toBeGreaterThan(0);
  const slices = await page.locator('#donutchart .nv-slice').count();
  const components = await page.evaluate(() => Object.values((window as any).selected_cc_count).filter((n) => (n as number) > 0).length);
  console.log(`donut chart: ${slices} slices for ${components} cellular components`);
  expect(slices).toBeGreaterThan(0);
  expect(slices).toBe(components);

  // attribute table after a click on the first node
  await openPanel(page, 'collapseTable');
  await page.evaluate(() => {
    const n = document.querySelector('.node')!;
    n.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true }));
  });
  await expect(page.locator('#rows tr')).not.toHaveCount(2);
  await expect(page.locator('#rows')).toContainText('Selected CC');
  await expect(page.locator('#rows')).toContainText('Extracted CCs');
  await expect(page.locator('#rows')).toContainText('InnateDB reference');

  // downloads
  await page.fill('#svg-name', 'my-network');
  const [svgDownload] = await Promise.all([page.waitForEvent('download'), page.click('#getpic')]);
  expect(svgDownload.suggestedFilename()).toBe('my-network.svg');
  const svgText = fs.readFileSync((await svgDownload.path())!, 'utf8');
  expect(svgText.startsWith('<svg')).toBe(true);
  expect(svgText).toMatch(/viewBox="0 0 850 850"/);
  expect(svgText).toMatch(/width="850"/);

  await page.selectOption('#figureformat', '.png');
  const [pngDownload] = await Promise.all([page.waitForEvent('download'), page.click('#getpic')]);
  expect(pngDownload.suggestedFilename()).toBe('my-network.png');
  const png = fs.readFileSync((await pngDownload.path())!);
  expect(png.subarray(1, 4).toString()).toBe('PNG');

  // CerebralWeb layout (menu item), at the end because it changes the page
  await page.click('a.dropdown-toggle:has-text("More options")');
  await page.click('#create_cerebral a');
  await expect(page.locator('#cerebralDiv')).toHaveCSS('visibility', 'visible');
  await expect.poll(() => page.locator('#cy2 canvas').count()).toBeGreaterThan(0);

  expect([...w.origins]).toEqual(['http://localhost:' + (process.env.PORT || 8080)]);
  expect(w.problems).toEqual([]);
});

test('responsive layout and screenshots', async ({ page }) => {
  const w = watch(page);
  for (const [width, height] of [
    [1280, 900],
    [390, 844],
  ]) {
    await page.setViewportSize({ width, height });
    await page.goto('/index.html');
    await loadExample(page);
    await page.waitForTimeout(3500); // let the layout settle for the picture
    const m = await page.evaluate(() => {
      const box = (s: string) => document.querySelector(s)!.getBoundingClientRect();
      return {
        scrollW: document.documentElement.scrollWidth,
        clientW: document.documentElement.clientWidth,
        diagram: box('#diagram'),
        left: box('.col-left'),
        content: box('#content'),
        right: box('.col-right'),
        banner: box('#legacy-banner'),
        notices: box('#notices'),
      };
    });
    expect(m.scrollW, `no horizontal overflow at ${width}px`).toBeLessThanOrEqual(m.clientW);
    expect(m.diagram.width).toBeLessThanOrEqual(850);
    expect(m.diagram.right).toBeLessThanOrEqual(m.clientW);
    expect(Math.abs(m.diagram.height - m.diagram.width)).toBeLessThan(2); // aspect ratio kept
    if (width < 992) {
      // panels stack: left column, then the diagram, then the attribute table
      expect(m.diagram.width).toBeGreaterThan(300);
      expect(m.content.top).toBeGreaterThanOrEqual(m.left.bottom - 1);
      expect(m.right.top).toBeGreaterThanOrEqual(m.content.bottom - 1);
      expect(Math.abs(m.left.left - m.content.left)).toBeLessThan(2);
    } else {
      expect(m.left.right).toBeLessThanOrEqual(m.content.left + 1);
      expect(m.content.right).toBeLessThanOrEqual(m.right.left + 1);
    }
    // the layout still works in 850 x 850 user units
    expect(await page.evaluate(() => (window as any).force.size())).toEqual([850, 850]);
    expect(await page.evaluate(() => document.getElementById('diagram')!.getAttribute('viewBox'))).toBe('0 0 850 850');

    // a click on a node still works on the scaled picture
    await openPanel(page, 'collapseTable');
    await page.evaluate(() => document.querySelector('.node')!.dispatchEvent(new MouseEvent('click', { bubbles: true })));
    await expect(page.locator('#rows')).toContainText('Selected CC');
    await openPanel(page, 'collapseSearch');
    await page.click('#clear');
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.screenshot({ path: path.join(root, 'test-results', 'screenshots', `example-${width}.png`), fullPage: width < 992 });
  }
  expect(w.problems).toEqual([]);
});

test('notices replace alert(): bad file, ?file= handling, edge bundling', async ({ page }) => {
  const w = watch(page);
  const base = 'http://localhost:' + (process.env.PORT || 8080);

  await page.goto('/index.html');
  await expect(page.locator('.group_component')).toHaveCount(20);

  // a file that is not XGMML
  await page.setInputFiles('#file', { name: 'bad.xgmml', mimeType: 'text/plain', buffer: Buffer.from('this is not xml') });
  const notice = page.locator('#notices .alert[role="status"]');
  await expect(notice).toHaveCount(1);
  await expect(notice).toContainText('could not be loaded');
  await expect(page.locator('#loading-container')).toBeHidden();
  await notice.locator('button.close').click();
  await expect(notice).toHaveCount(0);

  // edge bundling needs the forces stopped: a notice, not an alert
  await loadExample(page);
  await openPanel(page, 'collapseBundling');
  await page.click('#create_edge_bundling');
  await expect(page.locator('#notices')).toContainText("Can't compute edge bundling while Force algorithm is running");
  await openPanel(page, 'collapseForces');
  await page.click('#stop_forces');
  await page.click('#create_edge_bundling');
  await expect(page.locator('#notices')).toContainText('Edge bundling computed', { timeout: 60_000 });
  await page.click('#show_edge_bundling');
  await expect(page.locator('#diagram .bundling').first()).toBeAttached();
  await page.click('#remove_bundling');

  // ?file= : same-origin relative path
  await page.goto('/index.html?file=examples/rig-i-like-receptor-innatedb.xgmml');
  await expect(page.locator('#Nnodes')).toHaveText(`${EXAMPLE_NODES} nodes`);
  await expect(page.locator('#Nedges')).toHaveText(`${EXAMPLE_EDGES} links`);
  await expect(page.locator('#loading-container')).toBeHidden();

  // ?file= : only http(s)
  for (const bad of ['javascript:alert(1)', 'data:text/xml,<graph/>', 'file:///etc/passwd']) {
    await page.goto('/index.html?file=' + encodeURIComponent(bad));
    await expect(page.locator('#notices')).toContainText('not valid');
    await expect(page.locator('#Nnodes')).toHaveText('0 nodes');
    await expect(page.locator('.group_component')).toHaveCount(20); // the empty diagram is still shown
  }

  // ?file= : an encoded address on another origin that allows CORS (mocked), query string kept
  let requested = '';
  await page.route('https://data.example.test/**', (route) => {
    requested = route.request().url();
    route.fulfill({ status: 200, contentType: 'text/xml', headers: { 'access-control-allow-origin': '*' }, body: exampleXml });
  });
  await page.goto('/index.html?file=' + encodeURIComponent('https://data.example.test/net.xgmml?a=b'));
  await expect(page.locator('#Nnodes')).toHaveText(`${EXAMPLE_NODES} nodes`);
  expect(requested.startsWith('https://data.example.test/net.xgmml?a=b')).toBe(true);

  // ?file= : unreachable address -> notice, empty diagram
  await page.route('https://down.example.test/**', (route) => route.fulfill({ status: 500, body: 'no' }));
  await page.goto('/index.html?file=' + encodeURIComponent('https://down.example.test/x.xgmml'));
  await expect(page.locator('#notices')).toContainText('Error loading the XGMML file');
  await expect(page.locator('.group_component')).toHaveCount(20);
  await expect(page.locator('#loading-container')).toBeHidden();

  // the only requests to other origins were the ones mocked above
  const foreign = [...w.origins].filter((o) => o !== base);
  expect(foreign.sort()).toEqual(['https://data.example.test', 'https://down.example.test']);
  // (the 500 above is the only HTTP error and the console error that goes with it)
  expect(w.problems.filter((p) => !/down\.example\.test|status of 500|Failed to load resource/.test(p))).toEqual([]);
});

const NO_CC_XGMML = `<?xml version="1.0" encoding="UTF-8"?>
<graph label="no cc" xmlns="http://www.cs.rpi.edu/XGMML" directed="1">
  <node label="AAA" id="1"><att type="string" name="UniProt" value="P00001"/></node>
  <node label="BBB" id="2"><att type="string" name="UniProt" value="P00002"/></node>
  <node label="CCC" id="3"><att type="string" name="UniProt" value="P00003"/></node>
  <edge label="e1" source="1" target="2"/>
  <edge label="e2" source="2" target="3"/>
  <edge label="e3" source="3" target="999"/>
</graph>`;

test('InnateDB localization: asynchronous, https, clear message when unavailable', async ({ page }) => {
  const w = watch(page);
  const calls: string[] = [];
  let mode: 'ok' | 'down' = 'down';
  await page.route('https://www.innatedb.com/**', (route) => {
    calls.push(route.request().url());
    if (mode === 'down') return route.fulfill({ status: 503, body: 'down', headers: { 'access-control-allow-origin': '*' } });
    return route.fulfill({
      status: 200,
      contentType: 'application/json',
      headers: { 'access-control-allow-origin': '*' },
      body: JSON.stringify({ P00001: 'Nucleus', P00002: 'Cytosol', P00003: 'Plasma Membrane' }),
    });
  });
  const file = { name: 'nocc.xgmml', mimeType: 'text/xml', buffer: Buffer.from(NO_CC_XGMML) };

  await page.goto('/index.html');
  await expect(page.locator('.group_component')).toHaveCount(20);

  // service down: the question panel, then a clear notice, nothing drawn
  await page.setInputFiles('#file', file);
  await expect(page.locator('#cc_question')).toBeVisible();
  await page.click('#submit_cc');
  await expect(page.locator('#notices')).toContainText('InnateDB localization service');
  await expect(page.locator('#notices')).toContainText('unavailable');
  await expect(page.locator('#notices')).toContainText('"Selected CC"');
  await expect(page.locator('#loading-container')).toBeHidden();
  expect(calls.length).toBe(1);
  expect(calls[0].startsWith('https://www.innatedb.com/cerebralLocalizationWS.do?xref=UniProt')).toBe(true);

  // service back: load the same file twice, press Submit once: one request (handlers do not pile up)
  mode = 'ok';
  calls.length = 0;
  await page.setInputFiles('#file', file);
  await expect(page.locator('#cc_question')).toBeVisible();
  await page.setInputFiles('#file', { name: 'nocc2.xgmml', mimeType: 'text/xml', buffer: Buffer.from(NO_CC_XGMML) });
  await expect(page.locator('#cc_question')).toBeVisible();
  await page.click('#submit_cc');
  await expect(page.locator('#Nnodes')).toHaveText('3 nodes');
  await expect(page.locator('#Nedges')).toHaveText('2 links'); // the edge to a missing node is dropped
  expect(calls.length).toBe(1);
  await expect(page.locator('#color option[value="Selected CC"]')).toHaveCount(1);
  await expect(page.locator('#change_label option[value="Selected CC"]')).toHaveCount(1);
  expect(await page.evaluate(() => (window as any).d3.selectAll('.node').data().map((d: any) => d.component))).toEqual(['nucleus', 'cytoplasm', 'plasma_membrane']) // "cytosol" is drawn in the cytoplasm;

  expect(w.problems.filter((p) => !/503|Failed to load resource/.test(p))).toEqual([]);
});

test('attribute values from a file are shown as text, not as HTML', async ({ page }) => {
  const w = watch(page);
  const evil = `<?xml version="1.0"?>
<graph label="x" xmlns="http://www.cs.rpi.edu/XGMML">
  <node label="N1" id="1"><att type="string" name="Selected CC" value="nucleus"/><att type="string" name="Note" value="&lt;img src=x onerror=&quot;window.__pwned=1&quot;&gt;"/></node>
  <node label="N2" id="2"><att type="string" name="Selected CC" value="cytoplasm"/><att type="string" name="Note" value="plain"/></node>
</graph>`;
  await page.goto('/index.html');
  await expect(page.locator('.group_component')).toHaveCount(20);
  await page.setInputFiles('#file', { name: 'evil.xgmml', mimeType: 'text/xml', buffer: Buffer.from(evil) });
  await expect(page.locator('#Nnodes')).toHaveText('2 nodes');
  await openPanel(page, 'collapseTable');
  await page.evaluate(() => document.querySelector('.node')!.dispatchEvent(new MouseEvent('click', { bubbles: true })));
  await expect(page.locator('#rows')).toContainText('<img src=x');
  expect(await page.evaluate(() => (window as any).__pwned)).toBeUndefined();
  expect(await page.locator('#rows img').count()).toBe(0);
  expect(w.problems).toEqual([]);
});
