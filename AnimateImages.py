import cv2
import os
import re


def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)



# Image folder we are drawing our .pngs from
image_folder = 'Chirikov_Standard_Maps'

# Name of generated movie file
video_name = 'Chirikov_Standard_Maps.avi'

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]

sort_nicely(images)

frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 

video = cv2.VideoWriter(video_name, fourcc,frameSize= (width,height), fps=10)

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
