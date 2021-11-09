const puppeteer = require('puppeteer');
const util = require('../util');
it('restore_view"', async () => {
    const browser = await puppeteer.launch({headless: true});
    const page = await browser.newPage();
    await page.setViewport({width: 1500, height: 1000});
    await page.goto('http://127.0.0.1:5000/#q={"dataset":"../test-data/pbmc3k_no_raw.h5ad","q":[{"id":"n_counts","type":"obs"},{"id":"ABCB1","type":"X"}]}');
    await page.waitForSelector('[data-testid="scatter-chart-three"]');
    await page.evaluate(() => {
        document.querySelector('[data-testid="chart-extra"]').style.display = 'none';
    });
    const element = await page.$('[data-testid="scatter-chart-three"] > canvas');
    await element.screenshot({path: 'restore_view.png'});
    await browser.close();
    await util.diffImages('restore_view.png', 'screenshots/restore_view.png', 0);
});



