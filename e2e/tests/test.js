const puppeteer = require('puppeteer');
const gm = require('gm')
it('test"', async () => {

    const browser = await puppeteer.launch();
    const page = await browser.newPage();
    await page.goto('https://example.com');
    await page.waitForSelector('h1')
    await page.screenshot({path: 'example.png'});
    await browser.close();
});