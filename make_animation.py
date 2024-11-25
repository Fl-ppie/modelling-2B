import glob
    
from PIL import Image

def compare(x, y):
    if x=='final':
        return y
    if y=='final':
        return x
    
    if int(x)>int(y):
        return x
    if int(y)>int(x):
        return y
        

def make_gif():
    frames = [Image.open(image) for image in sorted(glob.glob("*.png"), key=compare)]
    
    frame_one=frames[0]
    frame_one.save("animation.gif",format='GIF', append_images=frames,
                   save_all=True, duration=100, loop=0)

if __name__ == "__main__":
    make_gif()