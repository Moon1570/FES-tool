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

  const CONFIG = {
    displaylogo: false,
    responsive: true,
    // Poster/slide export straight from the chart, since matplotlib is gone.
    toImageButtonOptions: { format: 'png', scale: 3, filename: 'fes-chart' },
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
  };

  function draw(el, traces, overrides) {
    if (!el) return;
    return Plotly.newPlot(el, traces, layout(overrides), CONFIG);
  }

  global.FES = { C, CYCLE_COLORS, FONT, layout, draw, CONFIG };
})(window);
