import time
import math
import random
import array
#from ulab import numpy
from cosmic import CosmicUnicorn, Channel
from picographics import PicoGraphics, DISPLAY_COSMIC_UNICORN as DISPLAY

'''
Population of 31 neurons randomly firing.
INSTALL ON PICO: Boot Pico in bootloader mode by holding down the BOOTSEL button while you plut into computer USB.
In Thonny, save as or drag and drop the file (rename it to to main.py) to the RPI_... drive that loaded up in bootloader mode.
THEN UNPLUG FROM THE COMPUTER!!!! I missed this part. If you plut back into a computer, then the main.py will be removed automatically and replaced by the default program.

You can adjust the cycling speed with A and B,
stripe width with C and D, hue with VOL + and -,
and the brightness with LUX + and -.
The sleep button stops the animation (can be started again with A or B).
'''
cu = CosmicUnicorn()
graphics = PicoGraphics(DISPLAY)

width = CosmicUnicorn.WIDTH
height = CosmicUnicorn.HEIGHT
notes = [random.randint(0, height) for _ in range(width)]


@micropython.native  # noqa: F821
def from_hsv(h, s, v):
    i = math.floor(h * 6.0)
    f = h * 6.0 - i
    v *= 255.0
    p = v * (1.0 - s)
    q = v * (1.0 - f * s)
    t = v * (1.0 - (1.0 - f) * s)

    i = int(i) % 6
    if i == 0:
        return int(v), int(t), int(p)
    if i == 1:
        return int(q), int(v), int(p)
    if i == 2:
        return int(p), int(v), int(t)
    if i == 3:
        return int(p), int(q), int(v)
    if i == 4:
        return int(t), int(p), int(v)
    if i == 5:
        return int(v), int(p), int(q)


@micropython.native  # noqa: F821
def draw():
    global hue_offset, phase
    phase_percent = phase / 15
    for x in range(width):
        colour = hue_map[int((x + (hue_offset * width)) % width)]
        for y in range(height):
            v = ((math.sin((x + y) / stripe_width + phase_percent) + 1.5) / 2.5)

            graphics.set_pen(graphics.create_pen(int(colour[0] * v), int(colour[1] * v), int(colour[2] * v)))
            graphics.pixel(x, y)

    cu.update(graphics)

boopety_beepety = cu.synth_channel(0)
boopety_beepety.configure(
    waveforms=Channel.SQUARE | Channel.SINE,
    attack=0.1,
    decay=0.2,
    sustain=0.0,
    release=0.5,
    volume=0.1
)
cu.play_synth()

hue_map = [from_hsv(x / width, 1.0, 1.0) for x in range(width)]
hue_offset = 0.0

play_sound = True

animate = True
stripe_width = 3.0
speed = 1.0
col_count = 0
next_col = 1
last_update_ms = 0.0
frate = array.array('f')
for iN in range(height):
    frate.append(random.uniform(.5,10))

cu.set_brightness(0.5)


def note_to_frequency(note_number):
    return int((2 ** ((note_number - 69.0) / 12)) * 440)

phase = 0
while True:

    if animate:
        phase += speed

    if cu.is_pressed(CosmicUnicorn.SWITCH_VOLUME_UP):
        play_sound = not play_sound

    if cu.is_pressed(CosmicUnicorn.SWITCH_VOLUME_DOWN):
        hue_offset -= 0.01
        hue_offset = 0.0 if hue_offset < 0.0 else hue_offset

    if cu.is_pressed(CosmicUnicorn.SWITCH_BRIGHTNESS_UP):
        cu.adjust_brightness(+0.01)

    if cu.is_pressed(CosmicUnicorn.SWITCH_BRIGHTNESS_DOWN):
        cu.adjust_brightness(-0.01)

    if cu.is_pressed(CosmicUnicorn.SWITCH_SLEEP):
        animate = False

    if cu.is_pressed(CosmicUnicorn.SWITCH_A):
        play_sound = not play_sound


    if cu.is_pressed(CosmicUnicorn.SWITCH_B):
        for iN in range(height):
            frate[iN] = frate[iN] + .01

    if cu.is_pressed(CosmicUnicorn.SWITCH_C):
        for iN in range(height):
            frate[iN] = frate[iN] - .01
            if frate[iN] < 0:
                frate[iN] = 0

    if cu.is_pressed(CosmicUnicorn.SWITCH_D):
        stripe_width -= 0.05
        stripe_width = 1.0 if stripe_width < 1.0 else stripe_width

    start = time.ticks_ms()
    
    #draw()
    # print(f'Hello, {os.getlogin()}! How are you?')
    if start-last_update_ms > 100:
        last_update_ms = start
        for iN in range(height):
            # Turn subsequent pixel off
            if random.uniform(.5,10) < frate[iN]/7:
                row_color = hue_map[int((iN + (hue_offset * height)) % height)]
                print(f'n{iN} f{frate[iN]} t{col_count} t{row_color}')
                graphics.set_pen(graphics.create_pen(int(row_color[0]), int(row_color[1]), int(row_color[2])))
                graphics.pixel(col_count, iN)
                if play_sound:
                    current_freq = note_to_frequency(iN)
                    boopety_beepety.frequency(current_freq+140)
                    boopety_beepety.trigger_attack()
                #graphics.set_pen(graphics.create_pen(int(255), int(0), int(0)))
                #graphics.pixel(next_col, iN)
                
        graphics.set_pen(graphics.create_pen(0, 0, 255))
        graphics.pixel(col_count, 31)
                
        cu.update(graphics)           
        
        print("total took: {} ms".format(time.ticks_ms() - start))
        col_count = col_count + 1
        next_col = col_count + 1
        
            
        
        if col_count == 32:
            col_count = 0
            # blank out screen
            graphics.set_pen(graphics.create_pen(int(0), int(0), int(0)))
            graphics.clear()
            
        if next_col == 32:
            next_col = 0
