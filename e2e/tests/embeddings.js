const puppeteer = require('puppeteer');
const gm = require('gm');

async function featureScreenshot(options) {
    const browser = await puppeteer.launch({headless: true});
    const page = await browser.newPage()
    await page.setViewport({width: 1500, height: 1000})
    await page.goto('http://localhost:3000/')

    await page.waitForSelector('[data-testid="' + options.input + '"]')
    await page.click('[data-testid="genes-input"]')

    await page.keyboard.type(options.name)
    await page.keyboard.press('Enter')
    await page.waitForSelector('[data-testid="scatter-chart-three"]')
    await page.waitForTimeout(100);
    await page.evaluate(() => {
        document.querySelector('[data-testid="chart-extra"]').remove();
    });
    const element = await page.$('[data-testid="scatter-chart-three"] > canvas');
    await element.screenshot({path: options.path})
    return {page, browser}
}

function diffImages(image1, image2, tolerance) {
    return new Promise((resolve, reject) => {
        gm.compare(image1, image2, tolerance, function (err, isEqual, equality, raw, path1, path2) {
            resolve(isEqual);
        })
    });
}

it('embeddings"', async () => {
    const {page, browser} = await featureScreenshot({name: 'CST3', input: 'genes-input', path: 'CST3.png'});
    const gallery = await page.$('[data-testid="gallery-image"]');
    await gallery.click();
    const element = await page.$('[data-testid="scatter-chart-three"] > canvas');
    await element.screenshot({path: 'louvain.png'})
    await page.waitForTimeout(100);
    await browser.close();
    await diffImages('CST3.png', 'screenshots/CST3.png', 0);
    // categories are drawn in random order
    await diffImages('louvain.png', 'screenshots/louvain.png', 0.01);
});



