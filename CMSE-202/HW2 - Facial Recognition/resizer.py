import PIL
from PIL import Image

import sys

img = Image.open(sys.argv[1])
img = img.resize((47, 62), PIL.Image.ANTIALIAS)
img.save(sys.argv[1])
