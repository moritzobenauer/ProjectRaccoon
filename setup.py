from setuptools import setup

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="raccoon",
    version="1.0.0",
    description="A Python program to generate 3D structures of linear polymers",
    author="Moritz L. Obenauer",
    url="https://github.com/moritzobenauer/ProjectRaccoon",
    packages=["raccoon"],
    package_dir={"raccoon": "src"},
    entry_points={"console_scripts": ["raccoon=raccoon.__main__:main"]},
    install_requires=required,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
