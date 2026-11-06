(() => {
  const menu = document.querySelector('.menu-toggle');
  const sidebar = document.querySelector('.sidebar');
  menu?.addEventListener('click', () => {
    const open = sidebar.classList.toggle('open');
    menu.setAttribute('aria-expanded', String(open));
  });
  sidebar?.addEventListener('click', event => {
    if (event.target.closest('a')) {
      sidebar.classList.remove('open');
      menu?.setAttribute('aria-expanded', 'false');
    }
  });
  const status = document.getElementById('copy-status');
  document.querySelectorAll('pre').forEach(pre => {
    const code = pre.querySelector('code');
    if (!code) return;
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'copy-button';
    button.textContent = 'Copy';
    button.setAttribute('aria-label', 'Copy code');
    button.addEventListener('click', async () => {
      try {
        if (!navigator.clipboard?.writeText) throw new Error('Clipboard unavailable');
        await navigator.clipboard.writeText(code.textContent);
        button.textContent = 'Copied';
        if (status) status.textContent = 'Code copied to clipboard.';
      } catch {
        const selection = window.getSelection();
        const range = document.createRange();
        range.selectNodeContents(code);
        selection.removeAllRanges();
        selection.addRange(range);
        button.textContent = 'Select';
        if (status) status.textContent = 'Code selected. Press Ctrl+C or Command+C to copy.';
      }
      setTimeout(() => { button.textContent = 'Copy'; }, 1800);
    });
    pre.append(button);
  });
  document.querySelectorAll('table').forEach(table => {
    if (table.parentElement.classList.contains('table-wrap')) return;
    const wrap = document.createElement('div');
    wrap.className = 'table-wrap';
    wrap.tabIndex = 0;
    wrap.setAttribute('role', 'region');
    wrap.setAttribute('aria-label', 'Scrollable reference table');
    table.before(wrap);
    wrap.append(table);
  });
  const navLinks = [...document.querySelectorAll('.sidebar a[href^="#"]')];
  const sections = [...document.querySelectorAll('section[id]')];
  if ('IntersectionObserver' in window) {
    const observer = new IntersectionObserver(entries => {
      const visible = entries.filter(entry => entry.isIntersecting).sort((a,b) => a.boundingClientRect.top-b.boundingClientRect.top);
      if (!visible.length) return;
      navLinks.forEach(link => {
        if (link.getAttribute('href') === '#' + visible[0].target.id) link.setAttribute('aria-current', 'location');
        else link.removeAttribute('aria-current');
      });
    }, {rootMargin:'-75px 0px -65% 0px'});
    sections.forEach(section => observer.observe(section));
  }
  const input = document.getElementById('doc-search');
  const results = document.getElementById('search-results');
  if (!input || !results) return;
  const searchable = [...document.querySelectorAll('[data-search]')].map(element => ({
    element,
    title: element.dataset.search,
    text: element.textContent.replace(/\s+/g, ' ').trim(),
    id: element.id
  }));
  function closeSearch() {
    results.hidden = true;
    input.setAttribute('aria-expanded', 'false');
  }
  input.addEventListener('input', () => {
    const query = input.value.trim().toLowerCase();
    results.replaceChildren();
    if (!query) { closeSearch(); return; }
    const words = query.split(/\s+/);
    const matches = searchable.filter(item => words.every(word => (item.title+' '+item.text).toLowerCase().includes(word)))
      .sort((a,b) => Number(b.title.toLowerCase().includes(query))-Number(a.title.toLowerCase().includes(query))).slice(0, 8);
    matches.forEach(item => {
      const link = document.createElement('a');
      link.href = '#' + item.id;
      link.textContent = item.title;
      const snippet = document.createElement('small');
      snippet.textContent = item.text.slice(0, 150) + (item.text.length > 150 ? '…' : '');
      link.append(snippet);
      link.addEventListener('click', closeSearch);
      results.append(link);
    });
    if (!matches.length) {
      const empty = document.createElement('p');
      empty.textContent = 'No matching section or parameter.';
      results.append(empty);
    }
    results.hidden = false;
    input.setAttribute('aria-expanded', 'true');
  });
  input.addEventListener('keydown', event => {
    if (event.key === 'Escape') closeSearch();
    if (event.key === 'ArrowDown') { results.querySelector('a')?.focus(); event.preventDefault(); }
    if (event.key === 'Enter') results.querySelector('a')?.click();
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape') { closeSearch(); sidebar?.classList.remove('open'); menu?.setAttribute('aria-expanded','false'); }
    if (event.key === '/' && !['INPUT','TEXTAREA'].includes(document.activeElement.tagName)) { event.preventDefault(); input.focus(); }
  });
  document.addEventListener('click', event => { if (!event.target.closest('.search')) closeSearch(); });
})();
