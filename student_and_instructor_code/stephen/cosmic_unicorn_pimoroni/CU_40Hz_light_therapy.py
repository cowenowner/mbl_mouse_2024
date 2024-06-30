import time
import math
import random
import array
from machine import Pin
#from ulab import numpy
from cosmic import CosmicUnicorn, Channel
from picographics import PicoGraphics, DISPLAY_COSMIC_UNICORN as DISPLAY

'''
Light stimulation for experiments and for therapy. In theory, and there is some date, 40Hz can triggery glymphatic amyloid beta removal.
https://www.cell.com/cell/fulltext/S0092-8674(19)30163-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419301631%3Fshowall%3Dtrue#secsectitle0080
mice receied stim for 1h/day for 7 days. It alternated 10s on and 10s off I believe.

INSTALL ON PICO: Boot Pico in bootloader mode by holding down the BOOTSEL button while you plug into computer USB. BOOTSEL is on the microcontroller, not the unicorn
In Thonny (and ONLY in THonny- copying in explorer won't work), save as and rename it to to main.py to the RPI_... 
THEN UNPLUG FROM THE COMPUTER!!!! I missed this part. If you plug back into a computer, then the main.py will be removed automatically and replaced by the default program.
for issues, see https://learn.pimoroni.com/article/getting-started-with-pico
'''
cu = CosmicUnicorn()
graphics = PicoGraphics(DISPLAY)

width = CosmicUnicorn.WIDTH
height = CosmicUnicorn.HEIGHT

#interval_sec = 0.0125
#interval_sec = 0.15
interval_sec = 0.0099
# Target is a delay of 25msec
volume_level = 0.05
#color_pen_background = graphics.create_pen(4, 253, 6)
color_pen_background = graphics.create_pen(0, 0, 0)
color_pen = graphics.create_pen(150, 10, 150)

boopety_beepety = cu.synth_channel(0)
boopety_beepety.configure(
    waveforms=Channel.SQUARE | Channel.SINE,
    attack=0.0001,
    decay=0.0001,
    sustain=0.1,
    release=0.1,
    volume=volume_level
)
tone1 = 2000 # The mice were at 10000Hz for 1ms
tone2 = 5100
play_sound = True



graphics.pixel(9, 29)
cu.set_brightness(0.5)
p0 = Pin(0, Pin.OUT)

boopety_beepety.play_tone(tone1, 1.0)



while True:
    start = time.ticks_ms()
    color_pen = graphics.create_pen(254,254,254)
    #color_pen = graphics.create_pen(random.randint(1,254), random.randint(1,254), random.randint(1,254))
    p0.value(0)
    if play_sound:
        boopety_beepety.volume(volume_level)
        time.sleep(0.001)
        boopety_beepety.volume(0.0)

        cu.play_synth()


    graphics.set_pen(color_pen_background)
    graphics.clear()
    cu.update(graphics)
    
    time.sleep(interval_sec)
    p0.value(1)
    #     if play_sound:
    #         #boopety_beepety.play_tone(tone2,.01)
    #         boopety_beepety.volume(0.0)
    #         cu.play_synth()
    
    
    graphics.set_pen(color_pen)
    graphics.clear()
    cu.update(graphics)
    
    time.sleep(interval_sec)
    
    print("took: {} ms".format(time.ticks_ms() - start))
    


