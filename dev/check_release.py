"""Check actual distribution members before release; standard library only.

Usage: python dev/check_release.py dist/*.whl dist/*.tar.gz
The development script itself is deliberately excluded from distributions.
"""
import argparse
from pathlib import Path, PurePosixPath
import re
import tarfile
import zipfile

RESOURCES = {'hastat/conf/GemmaPipelineConfig.toml', 'hastat/conf/bar.toml', 'hastat/conf/bar.mplstyle'}
METHODS = {f'{name}_METHODS.md' for name in ['COMPARE','FST','GENE_PLOT','HAPLOTYPE','LD','NETWORK','PLOT','POPULATION','STAT']}
SOURCE_ROOT = {'pyproject.toml','setup.py','setup.cfg','MANIFEST.in','README.md','CHANGELOG.md','LICENSE','PKG-INFO'} | METHODS
SOURCE_META = {'PKG-INFO','SOURCES.txt','dependency_links.txt','entry_points.txt','requires.txt','top_level.txt'}
WHEEL_META = {'METADATA','WHEEL','RECORD','entry_points.txt','top_level.txt','licenses/LICENSE'}
PACKAGE_DIRS = {'conf','dataio','func','log','utils','viz'}


def documentation_file(name):
    if name in {'docs/EXAMPLES.md','docs/gallery.html'}:
        return True
    if name.startswith('docs/images/'):
        leaf = name[len('docs/images/'):]
        return leaf in {case+'.'+ext for case in ['bar','pie','box','network','gene_tracks'] for ext in ['png','pdf','svg']}
    return name in {'docs/examples/'+leaf for leaf in ['gene.table.csv','gene.group.csv','populations.csv','phenotypes.csv',
        'genes.gff3','hap_colors.csv','annotations.csv','gene.pi.csv','gene.fst.csv','gene.network.csv',
        'gene_tracks.toml','commands.sh']}


def runtime_file(name):
    if name in RESOURCES:
        return True
    parts = PurePosixPath(name).parts
    if not parts or parts[0] != 'hastat' or not name.endswith('.py'):
        return False
    leaf = parts[-1]
    if leaf.startswith(('test_', 'validate')) or leaf.endswith('_test.py') or leaf == 'conftest.py':
        return False
    if len(parts) == 2:
        return leaf in {'__init__.py','hastat.py','manage.py'}
    return len(parts) == 3 and parts[1] in PACKAGE_DIRS


def audit(path):
    path = Path(path)
    if path.name.endswith('.whl'):
        match = re.fullmatch(r'hastat-([^-]+)-[^-]+-[^-]+-[^-]+\.whl',path.name)
        if not match:
            raise ValueError(f'Unexpected wheel name: {path.name}')
        version = match[1]
        meta = f'hastat-{version}.dist-info/'
        with zipfile.ZipFile(path) as archive:
            members = [item for item in archive.infolist() if not item.is_dir()]
            names = [item.filename for item in members]
            if any((item.external_attr >> 16) & 0o170000 == 0o120000 for item in members):
                raise ValueError('Symlink in wheel')
            payload = {name:archive.read(name) for name in names}
        allowed = lambda name: runtime_file(name) or name.startswith(meta) and name[len(meta):] in WHEEL_META
        required = RESOURCES | {'hastat/__init__.py','hastat/hastat.py'} | {meta+n for n in ['METADATA','WHEEL','RECORD','entry_points.txt','licenses/LICENSE']}
        metadata_name = meta+'METADATA'
        entry_name = meta+'entry_points.txt'
    elif path.name.endswith('.tar.gz'):
        match = re.fullmatch(r'hastat-(.+)\.tar\.gz',path.name)
        if not match:
            raise ValueError(f'Unexpected source archive name: {path.name}')
        version = match[1]
        root = f'hastat-{version}/'
        with tarfile.open(path) as archive:
            members = archive.getmembers()
            if any(not (item.isfile() or item.isdir()) for item in members):
                raise ValueError('Non-regular member in source archive')
            if any(not (item.name == root[:-1] or item.name.startswith(root)) for item in members):
                raise ValueError('Unexpected source archive root')
            names = [item.name[len(root):] for item in members if item.isfile()]
            payload = {item.name[len(root):]:archive.extractfile(item).read() for item in members if item.isfile()}
        meta = 'hastat.egg-info/'
        allowed = lambda name: runtime_file(name) or documentation_file(name) or name in SOURCE_ROOT or name.startswith(meta) and name[len(meta):] in SOURCE_META
        required = RESOURCES | (SOURCE_ROOT-{'setup.cfg'}) | {'hastat/__init__.py','hastat/hastat.py',meta+'entry_points.txt'}
        metadata_name = 'PKG-INFO'
        entry_name = meta+'entry_points.txt'
        # A stale egg-info file list can retain development files even when
        # the corresponding files were removed from the archive.
        for name in payload.get(meta+'SOURCES.txt',b'').decode().splitlines():
            if not allowed(name):
                raise ValueError(f'Non-release path in SOURCES.txt: {name}')
    else:
        raise ValueError(f'Unsupported distribution: {path.name}')
    if len(set(names)) != len(names):
        raise ValueError('Duplicate archive member')
    for name in names:
        parts = PurePosixPath(name).parts
        if name.startswith('/') or '\\' in name or '..' in parts or not allowed(name):
            raise ValueError(f'Non-release member: {name}')
    missing = required-set(names)
    if missing:
        raise ValueError(f'Missing required release files: {sorted(missing)}')
    metadata = payload[metadata_name].decode()
    if f'Version: {version}\n' not in metadata or 'License-Expression: BSD-3-Clause' not in metadata:
        raise ValueError('Version or license metadata mismatch')
    if 'hastat = hastat.hastat:main' not in payload[entry_name].decode():
        raise ValueError('Missing CLI entry point')
    return dict(version=version,files=sorted(names),runtime={name:body for name,body in payload.items() if runtime_file(name)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('distributions',nargs='+',type=Path)
    args = parser.parse_args()
    for path in args.distributions:
        result = audit(path)
        print(f'{path.name}: PASS ({len(result["files"])} files; no development/test payload)')
