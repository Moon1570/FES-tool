/* Shared Plotly configuration so every chart in the app reads as one system. */
(function (global) {
  'use strict';

  const C = {
    accent:  '#1f6feb',
    ok:      '#1a7f52',
    warn:    '#b26a00',
    danger:  '#c0392b',
    grey:    '#7a8796',
    grid:    '#e4e8ee',
    text:    '#11181f',
    text2:   '#4a5765',
  };

  // Categorical ramp for per-cycle series: cool -> warm as treatment progresses.
  const CYCLE_COLORS = [
    '#1f6feb', '#2b8ae0', '#22a0c8', '#1fae9e', '#3fae63',
    '#8aa63a', '#c09428', '#cf6f24', '#c0392b', '#a3306e',
    '#7b3fa8', '#4b47b5', '#2f6bd0', '#1f8fb0', '#26a37f',
    '#6fae42', '#b39a26', '#c8752a',
  ];

  const FONT = {
    family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif',
    size: 12,
    color: C.text2,
  };

  function layout(overrides) {
    const base = {
      font: FONT,
      margin: { l: 58, r: 18, t: 14, b: 44 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      hovermode: 'x unified',
      hoverlabel: { bgcolor: '#fff', bordercolor: C.grid, font: { size: 12 } },
      xaxis: { gridcolor: C.grid, zeroline: false, linecolor: C.grid, ticks: 'outside',
               tickcolor: C.grid, ticklen: 4 },
      yaxis: { gridcolor: C.grid, zeroline: false, linecolor: C.grid, ticks: 'outside',
               tickcolor: C.grid, ticklen: 4 },
      legend: { orientation: 'h', y: -0.2, x: 0, font: { size: 11 } },
      showlegend: false,
    };
    return Object.assign({}, base, overrides || {});
  }

  // Charts on the page stay simple: no toolbar, which would otherwise sit under the
  // full-screen button. The full-screen view is where people zoom, pan and export.
  const CONFIG = {
    displaylogo: false,
    responsive: true,
    displayModeBar: false,
  };

  const FULL_CONFIG = {
    displaylogo: false,
    responsive: true,
    displayModeBar: true,
    // Poster/slide export straight from the chart, since matplotlib is gone.
    toImageButtonOptions: { format: 'png', scale: 3, filename: 'fes-chart' },
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
  };

  const ICON_EXPAND = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" ' +
    'stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" ' +
    'aria-hidden="true"><path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/></svg>';
  const ICON_CLOSE = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" ' +
    'stroke="currentColor" stroke-width="2.4" stroke-linecap="round" aria-hidden="true">' +
    '<path d="M6 6l12 12M18 6L6 18"/></svg>';

  // ---------------------------------------------------------------- full screen
  // A page overlay rather than the browser Fullscreen API, which iPhones do not
  // support for anything but video. One overlay is shared by every chart.
  let overlay = null, fsChart = null, fsTitle = null, fsSub = null, returnFocus = null;

  function buildOverlay() {
    overlay = document.createElement('div');
    overlay.className = 'fs-overlay';
    overlay.hidden = true;
    overlay.setAttribute('role', 'dialog');
    overlay.setAttribute('aria-modal', 'true');
    overlay.innerHTML =
      '<div class="fs-head">' +
      '  <div class="fs-heading"><div class="fs-title"></div><div class="fs-sub"></div></div>' +
      '  <button type="button" class="fs-close" aria-label="Close full screen" title="Close (Esc)">' +
           ICON_CLOSE + '</button>' +
      '</div>' +
      '<div class="fs-chart"></div>' +
      '<div class="fs-hint">Drag across the chart to zoom in &middot; double-click to zoom out ' +
      '&middot; the toolbar pans and saves a picture</div>';
    document.body.appendChild(overlay);
    fsChart = overlay.querySelector('.fs-chart');
    fsTitle = overlay.querySelector('.fs-title');
    fsSub = overlay.querySelector('.fs-sub');
    overlay.querySelector('.fs-close').addEventListener('click', closeFullscreen);
    document.addEventListener('keydown', e => {
      if (e.key === 'Escape' && overlay && !overlay.hidden) closeFullscreen();
    });
  }

  function cardMeta(el) {
    const card = el.closest('.card');
    const h = card && card.querySelector(':scope > header h2');
    const hint = card && card.querySelector(':scope > header .hint');
    return { title: h ? h.textContent.trim() : 'Chart', subtitle: hint ? hint.textContent.trim() : '' };
  }

  function openFullscreen(el, opts, btn) {
    if (!el || !el.data || !el.data.length) return;
    opts = opts || {};
    if (!overlay) buildOverlay();
    const meta = (opts.meta && opts.meta()) || cardMeta(el);
    fsTitle.textContent = meta.title || 'Chart';
    fsSub.textContent = meta.subtitle || '';
    overlay.setAttribute('aria-label', meta.title || 'Chart');

    // Same traces and layout as the chart on the page, re-laid out for a big screen.
    const data = JSON.parse(JSON.stringify(el.data));
    const lay = JSON.parse(JSON.stringify(el.layout));
    delete lay.height;
    delete lay.width;
    lay.autosize = true;
    lay.margin = { l: 70, r: 30, t: 16, b: 60 };
    lay.font = Object.assign({}, FONT, { size: 14 });
    if (data.filter(t => t.name).length > 1) {
      lay.showlegend = true;
      // Legend on top on wide screens. On phones it wraps into several rows and would
      // run under Plotly's toolbar (always top-right), so it goes below the chart.
      lay.legend = window.matchMedia('(max-width: 700px)').matches
        ? { orientation: 'h', x: 0, y: -0.14, xanchor: 'left', yanchor: 'top', font: { size: 12 } }
        : { orientation: 'h', x: 0, y: 1.02, xanchor: 'left', yanchor: 'bottom', font: { size: 13 } };
    }
    if (opts.unifiedHover) lay.hovermode = 'x unified';

    returnFocus = btn || null;
    overlay.hidden = false;
    document.documentElement.classList.add('fs-open');
    Plotly.newPlot(fsChart, data, lay, FULL_CONFIG).then(() => {
      overlay.querySelector('.fs-close').focus();
    });
  }

  function closeFullscreen() {
    if (!overlay || overlay.hidden) return;
    Plotly.purge(fsChart);
    overlay.hidden = true;
    document.documentElement.classList.remove('fs-open');
    if (returnFocus) returnFocus.focus();
  }

  // Add the full-screen button for a chart. By default it goes at the right of the
  // chart's card header; pass opts.mount to place it elsewhere (the organ grid puts
  // it beside each organ's status badge so it never covers the data).
  function attachFullscreen(el, opts) {
    if (!el) return null;
    opts = opts || {};
    let btn = el._fesFsBtn;
    if (!btn) {
      btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'fs-btn';
      btn.title = 'Full screen';
      btn.setAttribute('aria-label', 'Show this chart full screen');
      btn.innerHTML = ICON_EXPAND;
      const card = el.closest('.card');
      const mount = opts.mount || (card && card.querySelector(':scope > header'));
      if (mount) {
        mount.appendChild(btn);
      } else {
        const wrap = document.createElement('div');
        wrap.className = 'fs-wrap';
        el.parentNode.insertBefore(wrap, el);
        wrap.appendChild(el);
        wrap.appendChild(btn);
        btn.classList.add('fs-float');
      }
      btn.addEventListener('click', () => openFullscreen(el, opts, btn));
      el._fesFsBtn = btn;
    }
    btn.hidden = !(el.data && el.data.length);
    return btn;
  }

  function draw(el, traces, overrides, opts) {
    if (!el) return;
    return Plotly.newPlot(el, traces, layout(overrides), CONFIG)
      .then(gd => { attachFullscreen(el, opts); return gd; });
  }

  global.FES = { C, CYCLE_COLORS, FONT, layout, draw, CONFIG, FULL_CONFIG,
                 attachFullscreen, openFullscreen, closeFullscreen };
})(window);
