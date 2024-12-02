import os
import glob

directory = 'images_newtonian_renamed/'

for path in glob.glob(directory+'*.png'):
    os.rename(path, path.replace(r'\0',r'\\').replace(r'\0',r'\\'))