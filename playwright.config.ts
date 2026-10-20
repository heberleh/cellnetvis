import { defineConfig } from '@playwright/test';

const port = Number(process.env.PORT) || 8080;

export default defineConfig({
  testDir: './tests',
  timeout: 60_000,
  expect: { timeout: 10_000 },
  fullyParallel: false,
  workers: 1,
  reporter: [['list']],
  use: {
    baseURL: `http://localhost:${port}`,
    acceptDownloads: true,
    launchOptions: {
      // Only needed when Playwright's own browser download is not installed (e.g. a sandbox with
      // PLAYWRIGHT_BROWSERS_PATH pointing at another revision): CHROMIUM_PATH=/path/to/chrome
      executablePath: process.env.CHROMIUM_PATH || undefined,
    },
  },
  webServer: {
    command: 'node scripts/serve.mjs',
    port,
    env: { PORT: String(port) },
    reuseExistingServer: true,
  },
});
