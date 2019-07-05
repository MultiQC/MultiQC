/////////////////////////////////////////////////
// Test if multiqc_report.html opens in a browser
/////////////////////////////////////////////////
// @ts-check

// With DEBUG=true, a browser window with devtools will open
var DEBUG = false;


var puppeteer = require('puppeteer');
var fs = require('fs');
// var beforeAll, describe, expect = require('jest');
jasmine.DEFAULT_TIMEOUT_INTERVAL = 60000;

/**
 * Headless browser handle
 * @type {import('puppeteer').Browser}
 */
var browser;

/**
 * Browser page handle
 * @type {import('puppeteer').Page}
 */
var page;

var browserErrors;
var silenceBrowserErrors;

beforeAll(async function () {
    browser = await puppeteer.launch({
        headless: !DEBUG,
        devtools: DEBUG,
        slowMo: DEBUG ? 1000 : 0,
        args: ['--no-sandbox'], // Needed on TravisCI
        defaultViewport: { width: 1280, height: 720 },
    });

    page = await browser.newPage();

    page.on('console', function (msg) {
        if (msg.type() === 'error') {
            browserErrors.push(msg.text());
        }
        if (!silenceBrowserErrors) console.log('page console:', msg.text());
    });
    page.on('pageerror', function (err) {
        browserErrors.push(err.toString());
        if (!silenceBrowserErrors) console.log('page error:', err.toString());
    });
    page.on('error', function (err) {
        browserErrors.push(err.toString());
        if (!silenceBrowserErrors) console.log('page error:', err.toString());
    });
});
beforeEach(function () {
    browserErrors = [];
    silenceBrowserErrors = false;
});


describe('Browser test', function () {
    describe('Content-Security-Policy', function () {
        beforeEach(async function () {
            // Clear page state before page.setContent(content)
            await page.goto('about:blank');
            await page.reload();
        });

        it('Opening multiqc_report.html with strict CSP should log violations', async function () {
            var originalContent = String(await fs.promises.readFile('./multiqc_report.html'));
            expect(originalContent).toContain('<head>');

            var content = originalContent.replace('<head>', `<head>
                <meta http-equiv="Content-Security-Policy" content="default-src 'self' 'report-sample';">
            `);

            expect(browserErrors).toEqual([]);
            silenceBrowserErrors = true;
            await page.setContent(content);
            silenceBrowserErrors = false;
            expect(browserErrors[0]).toContain('violates the following Content Security Policy');
        });

        it('CSP_HEADER.txt should contain up-to-date CSP hashes needed for multiqc_report.html', async function () {
            var originalContent = String(await fs.promises.readFile('./multiqc_report.html'));
            var cspHeader = String(await fs.promises.readFile('./../CSP_HEADER.txt'));
            var cleanCspHeader = cspHeader.replace(/#.*/g, '').replace(/[\n\r]/g, ' ');

            var content = originalContent.replace('<head>', `<head>
                <meta http-equiv="Content-Security-Policy" content="${cleanCspHeader}">
            `);

            expect(browserErrors).toEqual([]);
            await page.setContent(content);
            expect(browserErrors).toEqual([]);
        });
    });
});


afterAll(async function () {
    await page.close();
    await browser.close();
});
