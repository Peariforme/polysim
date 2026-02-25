/**
 * Shared navigation bar for benchmark pages.
 * Auto-injects a tab bar at the top of the page.
 * Used by both index.html (historical) and compare.html (comparison).
 */
(function () {
    'use strict';

    var path = window.location.pathname;
    var isCompare = path.indexOf('compare') !== -1;

    var nav = document.createElement('nav');
    nav.id = 'bench-nav';
    nav.innerHTML =
        '<div class="bench-nav-inner">' +
            '<span class="bench-nav-title">polysim</span>' +
            '<div class="bench-nav-tabs">' +
                '<a href="./index.html" class="bench-tab' + (isCompare ? '' : ' active') + '">' +
                    '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>' +
                    ' Historical Trends' +
                '</a>' +
                '<a href="./compare.html" class="bench-tab' + (isCompare ? ' active' : '') + '">' +
                    '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="3" y="12" width="4" height="9"/><rect x="10" y="5" width="4" height="16"/><rect x="17" y="9" width="4" height="12"/></svg>' +
                    ' Comparison' +
                '</a>' +
            '</div>' +
            '<a href="https://github.com/Peariforme/polysim" class="bench-nav-repo" title="GitHub">' +
                '<svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>' +
            '</a>' +
        '</div>';

    // Inject styles
    var style = document.createElement('style');
    style.textContent =
        '#bench-nav {' +
            'font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;' +
            'background: #1e293b;' +
            'border-bottom: 1px solid #334155;' +
            'padding: 0 1.5rem;' +
            'margin: -2rem -2rem 2rem -2rem;' +
            'position: sticky;' +
            'top: 0;' +
            'z-index: 1000;' +
        '}' +
        '.bench-nav-inner {' +
            'max-width: 1200px;' +
            'margin: 0 auto;' +
            'display: flex;' +
            'align-items: center;' +
            'gap: 1.5rem;' +
            'height: 3rem;' +
        '}' +
        '.bench-nav-title {' +
            'color: #e2e8f0;' +
            'font-weight: 700;' +
            'font-size: 0.9rem;' +
            'white-space: nowrap;' +
        '}' +
        '.bench-nav-tabs {' +
            'display: flex;' +
            'gap: 0.25rem;' +
            'flex: 1;' +
        '}' +
        '.bench-tab {' +
            'display: inline-flex;' +
            'align-items: center;' +
            'gap: 0.375rem;' +
            'padding: 0.5rem 0.875rem;' +
            'color: #94a3b8;' +
            'text-decoration: none;' +
            'font-size: 0.85rem;' +
            'font-weight: 500;' +
            'border-bottom: 2px solid transparent;' +
            'transition: color 0.15s, border-color 0.15s;' +
            'margin-bottom: -1px;' +
        '}' +
        '.bench-tab:hover {' +
            'color: #e2e8f0;' +
        '}' +
        '.bench-tab.active {' +
            'color: #60a5fa;' +
            'border-bottom-color: #3b82f6;' +
        '}' +
        '.bench-nav-repo {' +
            'color: #94a3b8;' +
            'display: inline-flex;' +
            'align-items: center;' +
            'transition: color 0.15s;' +
        '}' +
        '.bench-nav-repo:hover {' +
            'color: #e2e8f0;' +
        '}' +
        '@media (max-width: 480px) {' +
            '.bench-nav-title { display: none; }' +
            '.bench-tab { font-size: 0.8rem; padding: 0.5rem 0.5rem; }' +
        '}';

    // Wait for body if needed, then inject
    function inject() {
        document.body.insertBefore(nav, document.body.firstChild);
        document.head.appendChild(style);
    }

    if (document.body) {
        inject();
    } else {
        document.addEventListener('DOMContentLoaded', inject);
    }
})();
