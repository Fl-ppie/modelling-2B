import glob    
from PIL import Image

def make_gif(load_path, save_path, save_name):
    runtime = 30 # seconds
    
    # Load the images
    frames = [Image.open(image) for image in glob.glob(load_path+'*.png')]
    
    # Duplicate last image a few times
    for i in range(50):
        frames.append(frames[-1])
    
    # Create GIF
    frame_one=frames[0]
    frame_one.save(save_path+save_name+".gif",format='GIF', append_images=frames,
                   save_all=True, duration=runtime/len(frames)*1000, loop=0)
    
    print(f"Created GIF {save_name}.gif")

if __name__=='__main__':
    # Create a GIF per force type (might take a while)
    make_gif('z0.0086/position_data/images_newtonian/','','animation_newton')
    make_gif('z0.0086/position_data/images_mond/','','animation_mond')