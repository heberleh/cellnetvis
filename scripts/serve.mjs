// Tiny static file server for local development: serves web/ (no dependencies).
//   npm run dev            -> http://localhost:8080
//   PORT=3000 npm run dev  -> another port
// It also applies web/_headers (the same file Cloudflare Pages uses), so what you see locally
// (including Content-Security-Policy violations in the console) is what the deployed site does.
import http from 'node:http';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', 'web');
const port = Number(process.env.PORT) || 8080;

const types = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.gif': 'image/gif',
  '.jpg': 'image/jpeg',
  '.json': 'application/json; charset=utf-8',
  '.xgmml': 'application/xml; charset=utf-8',
  '.xml': 'application/xml; charset=utf-8',
  '.gpl': 'text/plain; charset=utf-8',
  '.woff': 'font/woff',
  '.woff2': 'font/woff2',
  '.ttf': 'font/ttf',
  '.eot': 'application/vnd.ms-fontobject',
};

// Minimal reader for the Cloudflare "_headers" format: a path pattern, then indented "Name: value" lines.
function readHeaders() {
  const rules = [];
  let current = null;
  let text = '';
  try {
    text = fs.readFileSync(path.join(root, '_headers'), 'utf8');
  } catch {
    return rules;
  }
  for (const line of text.split(/\r?\n/)) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    if (/^\s/.test(line)) {
      const i = line.indexOf(':');
      if (current && i > 0) current.headers[line.slice(0, i).trim()] = line.slice(i + 1).trim();
    } else {
      current = { pattern: line.trim(), headers: {} };
      rules.push(current);
    }
  }
  return rules;
}

function matches(pattern, urlPath) {
  const re = new RegExp('^' + pattern.replace(/[.+?^${}()|[\]\\]/g, '\\$&').replace(/\*/g, '.*') + '$');
  return re.test(urlPath);
}

const server = http.createServer((req, res) => {
  const url = new URL(req.url, 'http://localhost');
  let rel;
  try {
    rel = decodeURIComponent(url.pathname);
  } catch {
    res.writeHead(400).end('Bad request');
    return;
  }
  if (rel.endsWith('/')) rel += 'index.html';
  const file = path.normalize(path.join(root, rel));
  if (file !== root && !file.startsWith(root + path.sep)) {
    res.writeHead(403).end('Forbidden');
    return;
  }
  fs.readFile(file, (err, data) => {
    if (err) {
      res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' }).end('Not found');
      return;
    }
    const headers = { 'Content-Type': types[path.extname(file).toLowerCase()] || 'application/octet-stream' };
    for (const rule of readHeaders()) {
      if (matches(rule.pattern, url.pathname)) Object.assign(headers, rule.headers);
    }
    res.writeHead(200, headers).end(req.method === 'HEAD' ? undefined : data);
  });
});

server.listen(port, () => {
  console.log(`CellNetVis: http://localhost:${port}/  (serving ${root})`);
});
