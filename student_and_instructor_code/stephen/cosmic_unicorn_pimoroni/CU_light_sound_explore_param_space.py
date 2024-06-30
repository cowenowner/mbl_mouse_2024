import time
import math
import random
import array
from machine import Pin
#from ulab import numpy
from cosmic import CosmicUnicorn, Channel
from picographics import PicoGraphics, DISPLAY_COSMIC_UNICORN as DISPLAY

'''
For experiments: explore the response of neurons to a wide parameters space of sound and light.

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
volume_level = 0.25
pulse_count = 0
n_pulses_per_trial = 100
inter_trial_delay_sec = 5
trial_count = 0


condition_duration_ms = 5000
baseline_duration_ms = 5000
n_trials_per_condition = 12
stim_freq_to_explore = [1, 4, 10, 20, 40, 65]
light_cond_to_explore = [0, 1] # cycle between ligth on or off
sound_cond_to_explore = [0, 1] # cycle between sound on or off
#color_cond_to_explore = [0, 1] # cycle between changing color and constant color

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
time_in_trial_ms = 0 # keeps track of the cumulative time in trial.


graphics.pixel(9, 29)
cu.set_brightness(0.5)
p0 = Pin(2, Pin.OUT)
p1 = Pin(3, Pin.OUT)

boopety_beepety.play_tone(tone1, 1.0)



while True:
    start = time.ticks_ms()
    trial_start = time.ticks_ms()
    color_pen = graphics.create_pen(254,254,254)
    #color_pen = graphics.create_pen(random.randint(1,254), random.randint(1,254), random.randint(1,254))
    p0.value(0)
    p1.value(0)
    if play_sound:
        boopety_beepety.volume(volume_level)
        time.sleep(0.001)
        boopety_beepety.volume(0.0)

        cu.play_synth()
        p1.value(1)


    graphics.set_pen(color_pen_background)
    graphics.clear()
    cu.update(graphics)
    
    time.sleep(interval_sec)
    p0.value(1)
    #     if play_sound:
    #         #boopety_beepety.play_tone(tone2,.01)
    #         boopety_beepety.volume(0.0)
    #         cu.play_synth()
    if pulse_count >= n_pulses_per_trial:
        pulse_count = 0
        time.sleep(inter_trial_delay_sec)
        
    graphics.set_pen(color_pen)
    graphics.clear()
    cu.update(graphics)
    pulse_count = pulse_count + 1
    time.sleep(interval_sec)
    
    print("took: {} ms".format(time.ticks_ms() - start))
    time_in_trial_ms = time.ticks_ms() - trial_start;
    


