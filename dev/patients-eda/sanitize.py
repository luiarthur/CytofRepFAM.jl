import os

def sanitize(fpath):
    fpath_clean = fpath.lower()
    if fpath != fpath_clean:
        print('Sanitizing file: {}'.format(fpath))
        os.rename(fpath, fpath_clean)

if __name__ == '__main__':
    datadir = 'data/raw'

    for f in os.listdir(datadir):
        sanitize('{}/{}'.format(datadir, f))

