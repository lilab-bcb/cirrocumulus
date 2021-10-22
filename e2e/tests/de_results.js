const puppeteer = require('puppeteer');

it('de_results"', async () => {
    const browser = await puppeteer.launch({headless: true});
    const page = await browser.newPage();
    await page.setViewport({width: 1500, height: 1000});
    await page.goto('http://127.0.0.1:5000#q={"dataset":"../test-data/pbmc3k_no_raw.h5ad", "jobId":"cirro-rank_genes_groups"}');
    await page.waitForSelector('[data-testid="results-tab"]');
    await page.click('[data-testid="results-tab"]');
    await page.waitForSelector('[data-testid="dot-plot-table"]');
    await browser.close();
});



