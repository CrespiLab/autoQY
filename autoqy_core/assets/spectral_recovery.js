// Keep a single pending snapshot and serialize writes to avoid stale overwrites.
(() => {
  let pending = null;
  let timer = null;
  let saving = false;
  let hadData = false;
  const status = (message) => {
    window.dash_clientside.set_props('autosave-message', {children: message});
  };
  const flush = async () => {
    if (saving || !pending) return;
    const snapshot = pending;
    pending = null;
    saving = true;
    try {
      const response = await fetch('/_autoqy_spectral_autosave', {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(snapshot), signal: AbortSignal.timeout(30000)
      });
      if (response.status === 409) {
        pending = null;
        status('Another session was loaded. Reload this page before making changes.');
        return;
      }
      if (!response.ok) throw new Error('Autosave failed');
      status(pending ? 'Saving changes…'
        : 'Autosaved locally at ' + new Date().toLocaleTimeString() + '.');
    } catch (_) {
      pending = pending || snapshot;
      status('Autosave unavailable. Retrying; keep AutoQY open.');
    } finally {
      saving = false;
      if (pending) {
        clearTimeout(timer);
        timer = setTimeout(flush, 5000);
      }
    }
  };
  window.autoqySpectralAutosave = (ids, args) => {
    const [dataset, nmr, source, clearClicks, isomerRequest, generation, ...rest] = args;
    if (!dataset && !nmr && !hadData && !clearClicks) {
      return window.dash_clientside.no_update;
    }
    hadData = hadData || !!dataset || !!nmr;
    const settings = {};
    ids.forEach((id, index) => { settings[id] = rest[index]; });
    const [concentrations, paths, visibility, names, visibilityIds, nameIds]
      = rest.slice(ids.length);
    const packed = dataset ? Object.assign({}, dataset) : null;
    if (packed) {
      const count = packed.labels.length;
      if (concentrations && concentrations.length === count) {
        packed.concentrations = concentrations;
      }
      if (paths && paths.length === count) packed.path_lengths = paths;
      packed.legend_names = (packed.legend_names || packed.labels).slice();
      packed.legend_visibility = (packed.legend_visibility || Array(count).fill(true)).slice();
      (nameIds || []).forEach((id, index) => {
        if (id.index < count && names[index]) packed.legend_names[id.index] = names[index];
      });
      (visibilityIds || []).forEach((id, index) => {
        if (id.index < count) packed.legend_visibility[id.index] = (visibility[index] || []).includes('on');
      });
    }
    pending = {version: 1, dataset: packed, nmr: nmr, source: source,
      settings: settings, isomer_request: isomerRequest, generation: generation};
    status('Saving changes…');
    clearTimeout(timer);
    timer = setTimeout(flush, 1000);
    return window.dash_clientside.no_update;
  };
})();
