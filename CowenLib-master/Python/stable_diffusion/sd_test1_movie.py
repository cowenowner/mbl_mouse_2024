# from https://github.com/nateraw/stable-diffusion-videos
# Stable Diffusion models at https://huggingface.co/join?next=%2FCompVis%2Fstable-diffusion-v-1-4-original cowenow 43F-ew34A_4  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
# There does not seem to be a way to specify a start image instead of just a text prompt. 
from stable_diffusion_videos import StableDiffusionWalkPipeline
import torch

pipeline = StableDiffusionWalkPipeline.from_pretrained(
    "CompVis/stable-diffusion-v1-4",
    torch_dtype=torch.float16,
    revision="fp16",
).to("cuda")

video_path = pipeline.walk(
    prompts=['a cat', 'a dog'],
    seeds=[42, 1337],
    num_interpolation_steps=3,
    height=512,  # use multiples of 64 if > 512. Multiples of 8 if < 512.
    width=512,   # use multiples of 64 if > 512. Multiples of 8 if < 512.
    output_dir='dreams',        # Where images/videos will be saved
    name='animals_test',        # Subdirectory of output_dir where images/videos will be saved
    guidance_scale=8.5,         # Higher adheres to prompt more, lower lets model take the wheel
    num_inference_steps=50,     # Number of diffusion steps per image generated. 50 is good default
)

# Alternative - i think this allows you to start with an image that is in the .dreams direcotry.
video_path = pipeline.walk(
    prompts=['a handsome man', 'a handsome lion'],
    seeds=[42, 43],
    num_interpolation_steps=35,
    height=512,  # use multiples of 64 if > 512. Multiples of 8 if < 512.
    width=512,   # use multiples of 64 if > 512. Multiples of 8 if < 512.
    output_dir='dreams',        # Where images/videos will be saved
    name='animals_test2',        # Subdirectory of output_dir where images/videos will be saved
    guidance_scale=1.5,         # Higher adheres to prompt more, lower lets model take the wheel
    resume=True,         # In theory, this lets you take over from an image in the output directory - put an image in teh direcotry.
    num_inference_steps=15,     # Number of diffusion steps per image generated. 50 is good default
)

# ALTERNATIVE: Run locally
from stable_diffusion_videos import StableDiffusionWalkPipeline, Interface
import torch

pipeline = StableDiffusionWalkPipeline.from_pretrained(
    "CompVis/stable-diffusion-v1-4",
    torch_dtype=torch.float16,
    revision="fp16",
).to("cuda")

interface = Interface(pipeline)
interface.launch()