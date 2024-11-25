import glob    
from PIL import Image

path = 'images/'

def make_gif(path):
    frames = [Image.open(image) for image in glob.glob(path+'*.png')]
    
    for i in range(50):
        frames.append(frames[-1])
    
    frame_one=frames[0]
    frame_one.save(path+"animation"+".gif",format='GIF', append_images=frames,
                   save_all=True, duration=100, loop=0)

if __name__ == "__main__":
    make_gif(path)