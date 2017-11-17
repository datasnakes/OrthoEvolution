from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.Tools.utils import csvtolist

file = 'organisms_and_taxonomy_ids.csv'
ids = csvtolist(file, 'Taxonomy ID')

ncbiftp = NcbiFTPClient(email='shutchins2@umc.edu')
path = '/ddn/home3/vallender/databases/window-masker-files/'
ncbiftp.getwindowmaskerfiles(taxonomy_ids=ids, download_path=path)
